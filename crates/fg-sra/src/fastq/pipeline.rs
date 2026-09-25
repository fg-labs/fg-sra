//! Converting batches of spots to output bytes on worker threads, and writing them in order.
//!
//! The main thread dispatches fixed-size batches of spot ids. For each batch it hands every
//! output's writer a one-shot receiver, in dispatch order, before handing the batch to the
//! worker pool. A worker reads, routes and formats the batch's spots, compresses each
//! output's bytes into complete BGZF blocks (or leaves them plain), and sends them through
//! the batch's one-shot senders. Each writer waits on its receivers in order, so output is in
//! spot order whichever worker finishes first.

use std::fmt;
use std::fs::File;
use std::io::{BufWriter, ErrorKind, Write};
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result, anyhow};
use bgzf::{BGZF_BLOCK_SIZE, CompressionLevel, Compressor};
use crossbeam_channel::{Receiver, Sender, bounded};

use super::counts::SpotCounts;
use super::format::RecordFormatter;
use super::reader::SpotSource;
use super::spot::{SpotOutcome, SpotRouter};
use crate::progress::ProgressLogger;

/// Spots per batch. It is fixed, whatever the thread count, because BGZF blocks end only at
/// batch ends and at every `BGZF_BLOCK_SIZE` bytes within a batch; that keeps output bytes
/// identical at any thread count. A multiple of 8,192, the rows per blob in archives loaded
/// from FASTQ, so that no blob is decoded by two workers.
pub const BATCH_SPOTS: i64 = 16_384;

/// Stack for each thread that calls libncbi-vdb, whose schema evaluation recurses deeply on
/// some archives; fasterq-dump gives its threads 16 MiB after overflowing smaller ones. Rust's
/// default is 2 MiB.
pub const VDB_THREAD_STACK_BYTES: usize = 16 * 1024 * 1024;

/// Size of each writer's output buffer.
const WRITE_BUFFER_BYTES: usize = 256 * 1024;

/// Batches queued for each writer, per worker: enough for every worker to have one batch
/// in hand and one queued. More only holds converted bytes in memory when an output's
/// reader is slower than the workers.
const BATCHES_PER_WORKER: usize = 2;

/// Where an output goes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputTarget {
    /// Standard output (`-`).
    Stdout,
    /// A file, FIFO or device.
    Path(PathBuf),
}

impl OutputTarget {
    /// The target for a command-line path: `-` means stdout.
    pub fn from_arg(path: &Path) -> Self {
        if path == Path::new("-") { Self::Stdout } else { Self::Path(path.to_path_buf()) }
    }

    /// A name for messages.
    pub fn describe(&self) -> String {
        match self {
            Self::Stdout => "stdout".to_owned(),
            Self::Path(path) => path.display().to_string(),
        }
    }
}

/// How an output's bytes are encoded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    /// Uncompressed text.
    Plain,
    /// BGZF blocks.
    Bgzf,
}

/// One output file (or stdout).
#[derive(Debug, Clone)]
pub struct Output {
    pub target: OutputTarget,
    pub encoding: Encoding,
}

/// Which output each kind of record goes to, as indices into the outputs.
#[derive(Debug, Clone, Default)]
pub struct OutputLayout {
    pub first_mates: Option<usize>,
    pub second_mates: Option<usize>,
    /// Both mates of each pair, first then second.
    pub interleaved: Option<usize>,
    pub unpaired: Option<usize>,
    /// Output of each technical read number, from 1.
    pub technical: Vec<usize>,
}

/// Everything about a conversion that workers need.
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    pub router: SpotRouter,
    pub formatter: RecordFormatter,
    pub layout: OutputLayout,
    /// Whether qualities are read, so spots must carry one per base.
    pub with_qualities: bool,
    /// BGZF compression level, for BGZF outputs.
    pub compression_level: CompressionLevel,
    /// Spots per batch: [`BATCH_SPOTS`], other than in tests.
    pub batch_spots: i64,
}

/// What a finished conversion counted.
#[derive(Debug, Default)]
pub struct PipelineSummary {
    /// Every spot converted, summed over the workers.
    pub counts: SpotCounts,
    /// Technical reads written, over every technical output.
    pub technical_reads_written: u64,
    /// An output's reader went away (a closed pipe), so conversion stopped early.
    pub stopped_early: bool,
}

/// A batch of spot ids, with one sender per output for the batch's encoded bytes.
struct Job {
    /// First spot id of the batch.
    first: i64,
    /// Last spot id of the batch (inclusive).
    last: i64,
    /// One sender per output, in output order.
    results: Vec<Sender<Result<Vec<u8>>>>,
}

impl Job {
    /// Tell each output's writer that this batch won't be converted.
    fn abandon(self) {
        for result in self.results {
            // A writer that has already exited no longer needs to hear.
            let _ = result.send(Err(Abandoned.into()));
        }
    }
}

/// Sent in place of a batch's bytes when the run is stopping; never the error reported.
#[derive(Debug)]
struct Abandoned;

impl fmt::Display for Abandoned {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "batch abandoned because the conversion is stopping")
    }
}

impl std::error::Error for Abandoned {}

/// One worker's counts.
#[derive(Default)]
struct WorkerTally {
    counts: SpotCounts,
    technical_reads_written: u64,
}

/// Convert the spots in `spots` to `outputs`, with one worker thread per source.
///
/// On failure, the first real error is returned: a worker's (reading, validating or
/// compressing) in preference to a writer's. If an output's reader goes away the run stops
/// early without error, and the summary says so. BGZF outputs get their end-of-file block
/// only when every spot was converted.
pub fn run<S: SpotSource + Send>(
    sources: Vec<S>,
    spots: RangeInclusive<i64>,
    outputs: &[Output],
    config: &PipelineConfig,
    progress: &ProgressLogger,
) -> Result<PipelineSummary> {
    let stop = AtomicBool::new(false);
    let pipe_closed = AtomicBool::new(false);
    let threads = sources.len().max(1);

    std::thread::scope(|scope| {
        let (job_tx, job_rx) = bounded::<Job>(threads);
        let mut order_txs = Vec::with_capacity(outputs.len());
        let mut writers = Vec::with_capacity(outputs.len());
        for output in outputs {
            let (order_tx, order_rx) = bounded(threads * BATCHES_PER_WORKER);
            order_txs.push(order_tx);
            let (stop, pipe_closed) = (&stop, &pipe_closed);
            writers.push(scope.spawn(move || write_output(&order_rx, output, stop, pipe_closed)));
        }
        let workers = sources
            .into_iter()
            .map(|source| {
                let (job_rx, stop) = (job_rx.clone(), &stop);
                std::thread::Builder::new()
                    .stack_size(VDB_THREAD_STACK_BYTES)
                    .spawn_scoped(scope, move || {
                        work(source, &job_rx, outputs, config, stop, progress)
                    })
            })
            .collect::<std::io::Result<Vec<_>>>()
            .context("failed to start a worker thread")?;
        drop(job_rx);

        dispatch(spots, config.batch_spots, &job_tx, &order_txs, &stop);
        drop(job_tx);
        drop(order_txs);

        let mut summary = PipelineSummary::default();
        let mut first_error = None;
        for worker in workers {
            match worker.join().unwrap_or_else(|_| Err(anyhow!("a worker thread panicked"))) {
                Ok(tally) => {
                    summary.counts.merge(&tally.counts);
                    summary.technical_reads_written += tally.technical_reads_written;
                }
                Err(e) => keep_first_real_error(&mut first_error, e),
            }
        }
        for writer in writers {
            if let Err(e) =
                writer.join().unwrap_or_else(|_| Err(anyhow!("a writer thread panicked")))
            {
                keep_first_real_error(&mut first_error, e);
            }
        }
        if let Some(e) = first_error {
            return Err(e);
        }
        summary.stopped_early = pipe_closed.load(Ordering::Relaxed);
        Ok(summary)
    })
}

/// Keep `error` as the run's error unless one is kept already or it is only [`Abandoned`].
fn keep_first_real_error(first: &mut Option<anyhow::Error>, error: anyhow::Error) {
    if first.is_none() && !error.is::<Abandoned>() {
        *first = Some(error);
    }
}

/// Hand out `spots` in batches until done or told to stop. Each writer is given its receiver
/// before the batch goes to a worker, so no batch's bytes can arrive unexpected.
fn dispatch(
    spots: RangeInclusive<i64>,
    batch_spots: i64,
    job_tx: &Sender<Job>,
    order_txs: &[Sender<Receiver<Result<Vec<u8>>>>],
    stop: &AtomicBool,
) {
    let (mut first, end) = spots.into_inner();
    while first <= end && !stop.load(Ordering::Relaxed) {
        let last = end.min(first + batch_spots.max(1) - 1);
        let mut results = Vec::with_capacity(order_txs.len());
        for order_tx in order_txs {
            let (result_tx, result_rx) = bounded(1);
            if order_tx.send(result_rx).is_err() {
                // The writer has failed; its error is reported when it is joined. The writers
                // already given this batch skip it rather than report it as never converted.
                stop.store(true, Ordering::Relaxed);
                Job { first, last, results }.abandon();
                return;
            }
            results.push(result_tx);
        }
        if let Err(unsent) = job_tx.send(Job { first, last, results }) {
            unsent.into_inner().abandon();
            return;
        }
        first = last + 1;
    }
}

/// Worker loop: convert each batch received, or abandon it once the run is stopping.
fn work<S: SpotSource>(
    mut source: S,
    jobs: &Receiver<Job>,
    outputs: &[Output],
    config: &PipelineConfig,
    stop: &AtomicBool,
    progress: &ProgressLogger,
) -> Result<WorkerTally> {
    let num_outputs = outputs.len();
    let mut compressors: Vec<_> =
        (0..num_outputs).map(|_| Compressor::new(config.compression_level)).collect();
    let mut buffers = vec![Vec::new(); num_outputs];
    let mut tally = WorkerTally::default();
    let mut failure = None;

    for job in jobs {
        if failure.is_some() || stop.load(Ordering::Relaxed) {
            job.abandon();
            continue;
        }
        buffers.iter_mut().for_each(Vec::clear);
        let converted = convert_batch(&mut source, job.first..=job.last, config, &mut buffers)
            .and_then(|batch_tally| {
                let encoded = encode(&mut compressors, &mut buffers, outputs)?;
                Ok((batch_tally, encoded))
            });
        match converted {
            Ok((batch_tally, encoded)) => {
                tally.counts.merge(&batch_tally.counts);
                tally.technical_reads_written += batch_tally.technical_reads_written;
                for (result, bytes) in job.results.into_iter().zip(encoded) {
                    // A writer that has already exited no longer needs the bytes.
                    let _ = result.send(Ok(bytes));
                }
                progress.record((job.last - job.first + 1) as u64);
            }
            Err(e) => {
                stop.store(true, Ordering::Relaxed);
                job.abandon();
                failure = Some(e);
            }
        }
    }
    failure.map_or(Ok(tally), Err)
}

/// Read, route and format the spots `ids` into one buffer per output.
fn convert_batch<S: SpotSource>(
    source: &mut S,
    ids: RangeInclusive<i64>,
    config: &PipelineConfig,
    buffers: &mut [Vec<u8>],
) -> Result<WorkerTally> {
    let layout = &config.layout;
    let formatter = &config.formatter;
    let mut tally = WorkerTally::default();
    for id in ids {
        let spot = source.read(id)?;
        spot.validate(config.with_qualities)?;
        let outcome = config.router.route(&spot);
        tally.counts.record(&spot, outcome);
        match outcome {
            SpotOutcome::Pair(first, second) => {
                if let Some(index) = layout.interleaved {
                    formatter.write(&mut buffers[index], &spot, first);
                    formatter.write(&mut buffers[index], &spot, second);
                } else if let (Some(first_index), Some(second_index)) =
                    (layout.first_mates, layout.second_mates)
                {
                    formatter.write(&mut buffers[first_index], &spot, first);
                    formatter.write(&mut buffers[second_index], &spot, second);
                }
            }
            SpotOutcome::Unpaired(read) => {
                if let Some(index) = layout.unpaired {
                    formatter.write(&mut buffers[index], &spot, read);
                }
            }
            SpotOutcome::Dropped(_) => continue,
        }
        // The router has checked a written spot has one technical read per technical output.
        for (read, &index) in spot.technical_reads().zip(&layout.technical) {
            formatter.write(&mut buffers[index], &spot, read);
            tally.technical_reads_written += 1;
        }
    }
    Ok(tally)
}

/// Encode each output's buffer: as complete BGZF blocks, or moved out as it is.
fn encode(
    compressors: &mut [Compressor],
    buffers: &mut [Vec<u8>],
    outputs: &[Output],
) -> Result<Vec<Vec<u8>>> {
    buffers
        .iter_mut()
        .zip(compressors)
        .zip(outputs)
        .map(|((buffer, compressor), output)| match output.encoding {
            Encoding::Plain => Ok(std::mem::take(buffer)),
            Encoding::Bgzf => {
                let mut compressed = Vec::with_capacity(buffer.len() / 3 + 1024);
                let mut block = Vec::with_capacity(BGZF_BLOCK_SIZE);
                for chunk in buffer.chunks(BGZF_BLOCK_SIZE) {
                    compressor
                        .compress(chunk, &mut block)
                        .map_err(|e| anyhow!("BGZF compression failed: {e}"))?;
                    compressed.extend_from_slice(&block);
                }
                Ok(compressed)
            }
        })
        .collect()
}

/// Writer loop for one output: write each batch's bytes in dispatch order.
///
/// A closed pipe (the output's reader went away) stops the whole run: this writer notes it,
/// then keeps draining its batches without writing them so no worker waits on it. The BGZF
/// end-of-file block is written only if the run wasn't stopped, so an incomplete output is
/// detectably truncated.
fn write_output(
    order: &Receiver<Receiver<Result<Vec<u8>>>>,
    output: &Output,
    stop: &AtomicBool,
    pipe_closed: &AtomicBool,
) -> Result<()> {
    let name = output.target.describe();
    let sink: Box<dyn Write> = match &output.target {
        OutputTarget::Stdout => Box::new(std::io::stdout().lock()),
        OutputTarget::Path(path) => Box::new(
            File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
        ),
    };
    let mut writer = BufWriter::with_capacity(WRITE_BUFFER_BYTES, sink);
    let mut closed = false;
    let note_closed = |closed: &mut bool| {
        *closed = true;
        pipe_closed.store(true, Ordering::Relaxed);
        stop.store(true, Ordering::Relaxed);
    };

    for batch in order {
        let bytes = match batch.recv() {
            Ok(Ok(bytes)) => bytes,
            Ok(Err(e)) if e.is::<Abandoned>() => continue,
            Ok(Err(e)) => return Err(e),
            Err(_) => return Err(anyhow!("a worker exited without converting a batch")),
        };
        if closed {
            continue;
        }
        if let Err(e) = writer.write_all(&bytes) {
            if e.kind() != ErrorKind::BrokenPipe {
                return Err(e).with_context(|| format!("failed to write {name}"));
            }
            note_closed(&mut closed);
        }
    }
    if closed || stop.load(Ordering::Relaxed) {
        return Ok(());
    }
    if output.encoding == Encoding::Bgzf {
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        if let Err(e) = writer.write_all(&eof) {
            if e.kind() != ErrorKind::BrokenPipe {
                return Err(e).with_context(|| format!("failed to write {name}"));
            }
            note_closed(&mut closed);
            return Ok(());
        }
    }
    match writer.flush() {
        Err(e) if e.kind() == ErrorKind::BrokenPipe => {
            note_closed(&mut closed);
            Ok(())
        }
        result => result.with_context(|| format!("failed to write {name}")),
    }
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::*;
    use crate::fastq::defline::Defline;
    use crate::fastq::spot::Spot;
    use crate::fastq::spot::tests::TestSpot;
    use crate::record::READ_TYPE_BIOLOGICAL as B;

    const T: u8 = 0;

    /// Serves test spots by id, 1-based; `fail_at` makes reading that spot an error.
    struct TestSource {
        spots: Vec<TestSpot>,
        fail_at: Option<i64>,
    }

    impl SpotSource for TestSource {
        fn read(&mut self, id: i64) -> Result<Spot<'_>> {
            if self.fail_at == Some(id) {
                anyhow::bail!("test read failure at spot {id}");
            }
            Ok(self.spots[(id - 1) as usize].view())
        }
    }

    /// A source whose every read first uses 8 MiB of stack, as libncbi-vdb can.
    struct DeepStackSource(TestSource);

    impl SpotSource for DeepStackSource {
        // A large stack frame is the point.
        #[allow(clippy::large_stack_arrays)]
        fn read(&mut self, id: i64) -> Result<Spot<'_>> {
            let stack = std::hint::black_box([1u8; 8 << 20]);
            assert_eq!(stack[stack.len() - 1], 1);
            self.0.read(id)
        }
    }

    /// `count` spots cycling through pairs, orphans of an empty first mate, and single
    /// reads, each with a leading technical read, and bases varying by spot.
    fn test_spots(count: usize) -> Vec<TestSpot> {
        (0..count)
            .map(|i| {
                let bases = ["ACGTACGTAC", "GGCCTTAAGGCCT", "TTTTGGGGCCCCAAAA"][i % 3];
                let (a, b) = bases.split_at(i % 5 + 2);
                let reads: Vec<(&str, u8, u8)> = match i % 3 {
                    0 => vec![("NN", T, 0), (a, B, 0), (b, B, 0)],
                    1 => vec![("NN", T, 0), ("", B, 0), (b, B, 0)],
                    _ => vec![("NN", T, 0), (a, B, 0)],
                };
                let mut spot = TestSpot::new(&reads);
                spot.id = i as i64 + 1;
                spot
            })
            .collect()
    }

    /// One source of `count` test spots per thread.
    fn sources(count: usize, threads: usize) -> Vec<TestSource> {
        (0..threads).map(|_| TestSource { spots: test_spots(count), fail_at: None }).collect()
    }

    /// Sources whose reads fail at spot `fail_at`.
    fn failing_sources(count: usize, threads: usize, fail_at: i64) -> Vec<TestSource> {
        let mut sources = sources(count, threads);
        for source in &mut sources {
            source.fail_at = Some(fail_at);
        }
        sources
    }

    fn config(layout: OutputLayout) -> PipelineConfig {
        PipelineConfig {
            router: SpotRouter {
                min_read_len: 0,
                kept_filters: [true; 4],
                write_pairs: layout.first_mates.is_some() || layout.interleaved.is_some(),
                write_unpaired: layout.unpaired.is_some(),
                technical_reads: (!layout.technical.is_empty()).then_some(layout.technical.len()),
            },
            formatter: RecordFormatter::new(Defline::parse("$ac.$si/$ri").unwrap(), "T", false),
            layout,
            with_qualities: true,
            compression_level: CompressionLevel::new(1).unwrap(),
            batch_spots: 4,
        }
    }

    /// A fresh directory for one test's outputs.
    fn scratch_dir(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("fg_sra_pipeline_{}_{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn output(dir: &Path, name: &str, encoding: Encoding) -> Output {
        Output { target: OutputTarget::Path(dir.join(name)), encoding }
    }

    fn paired_layout() -> OutputLayout {
        OutputLayout {
            first_mates: Some(0),
            second_mates: Some(1),
            unpaired: Some(2),
            ..OutputLayout::default()
        }
    }

    fn run_to<S: SpotSource + Send>(
        sources: Vec<S>,
        spots: RangeInclusive<i64>,
        outputs: &[Output],
        config: &PipelineConfig,
    ) -> Result<PipelineSummary> {
        run(sources, spots, outputs, config, &ProgressLogger::new(0, 0))
    }

    fn read_text(path: &Path) -> String {
        std::fs::read_to_string(path).unwrap()
    }

    fn names(fastq: &str) -> Vec<&str> {
        fastq.lines().step_by(4).collect()
    }

    #[test]
    fn mates_go_to_their_outputs_in_spot_order() {
        let dir = scratch_dir("mates");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        run_to(sources(9, 3), 1..=9, &outputs, &config(paired_layout())).unwrap();
        assert_eq!(names(&read_text(&dir.join("r1.fq"))), ["@T.1/1", "@T.4/1", "@T.7/1"]);
        assert_eq!(names(&read_text(&dir.join("r2.fq"))), ["@T.1/2", "@T.4/2", "@T.7/2"]);
        let unpaired = read_text(&dir.join("u.fq"));
        assert_eq!(names(&unpaired), ["@T.2/2", "@T.3/1", "@T.5/2", "@T.6/1", "@T.8/2", "@T.9/1"]);
    }

    #[test]
    fn workers_have_room_for_deep_library_stacks() {
        let dir = scratch_dir("deep-stack");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let sources: Vec<_> = sources(3, 2).into_iter().map(DeepStackSource).collect();
        run_to(sources, 1..=3, &outputs, &config(paired_layout())).unwrap();
        assert_eq!(names(&read_text(&dir.join("r1.fq"))), ["@T.1/1"]);
    }

    #[test]
    fn only_the_requested_spots_are_converted() {
        let dir = scratch_dir("range");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let summary = run_to(sources(9, 2), 4..=6, &outputs, &config(paired_layout())).unwrap();
        assert_eq!(names(&read_text(&dir.join("r1.fq"))), ["@T.4/1"]);
        assert_eq!(summary.counts.spots, 3);
    }

    #[test]
    fn interleaved_output_holds_each_pairs_mates_in_turn() {
        let dir = scratch_dir("interleaved");
        let outputs = [output(&dir, "pairs.fq", Encoding::Plain)];
        let layout = OutputLayout { interleaved: Some(0), ..OutputLayout::default() };
        run_to(sources(7, 2), 1..=7, &outputs, &config(layout)).unwrap();
        let names = names(&read_text(&dir.join("pairs.fq"))).join(" ");
        assert_eq!(names, "@T.1/1 @T.1/2 @T.4/1 @T.4/2 @T.7/1 @T.7/2");
    }

    #[test]
    fn technical_reads_line_up_with_the_spots_written() {
        let dir = scratch_dir("technical");
        let outputs = ["u.fq", "t1.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let layout =
            OutputLayout { unpaired: Some(0), technical: vec![1], ..OutputLayout::default() };
        let summary = run_to(sources(6, 2), 1..=6, &outputs, &config(layout)).unwrap();
        let unpaired: Vec<String> = names(&read_text(&dir.join("u.fq")))
            .iter()
            .map(|n| n.split('/').next().unwrap().to_owned())
            .collect();
        let technical: Vec<String> = names(&read_text(&dir.join("t1.fq")))
            .iter()
            .map(|n| n.split('/').next().unwrap().to_owned())
            .collect();
        assert_eq!(unpaired, technical);
        assert_eq!(summary.technical_reads_written, technical.len() as u64);
    }

    #[test]
    fn counts_cover_every_spot_converted() {
        let dir = scratch_dir("counts");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let summary = run_to(sources(10, 4), 1..=10, &outputs, &config(paired_layout())).unwrap();
        assert_eq!(summary.counts.spots, 10);
        assert_eq!(summary.counts.pairs_written, 4);
        assert_eq!(summary.counts.unpaired_written, 6);
        assert!(!summary.stopped_early);
    }

    #[test]
    fn output_bytes_are_the_same_at_any_thread_count() {
        let convert = |threads: usize| {
            let dir = scratch_dir(&format!("threads{threads}"));
            let outputs = [
                output(&dir, "r1.fq.gz", Encoding::Bgzf),
                output(&dir, "r2.fq", Encoding::Plain),
                output(&dir, "u.fq.gz", Encoding::Bgzf),
            ];
            run_to(sources(2_000, threads), 1..=2_000, &outputs, &config(paired_layout())).unwrap();
            ["r1.fq.gz", "r2.fq", "u.fq.gz"].map(|n| std::fs::read(dir.join(n)).unwrap())
        };
        let single = convert(1);
        assert_eq!(single, convert(3));
        assert_eq!(single, convert(8));
    }

    #[test]
    fn bgzf_output_decompresses_to_the_plain_output() {
        let dir = scratch_dir("bgzf");
        let outputs = [
            output(&dir, "r1.fq.gz", Encoding::Bgzf),
            output(&dir, "r2.fq", Encoding::Plain),
            output(&dir, "u.fq", Encoding::Plain),
        ];
        let layout = paired_layout();
        run_to(sources(3_000, 4), 1..=3_000, &outputs, &config(layout.clone())).unwrap();
        let plain_outputs = ["p1.fq", "p2.fq", "pu.fq"].map(|n| output(&dir, n, Encoding::Plain));
        run_to(sources(3_000, 4), 1..=3_000, &plain_outputs, &config(layout)).unwrap();

        let mut decompressed = String::new();
        flate2::read::MultiGzDecoder::new(File::open(dir.join("r1.fq.gz")).unwrap())
            .read_to_string(&mut decompressed)
            .unwrap();
        assert_eq!(decompressed, read_text(&dir.join("p1.fq")));
    }

    #[test]
    fn bgzf_output_ends_with_the_end_of_file_block() {
        let dir = scratch_dir("eof");
        let outputs = [output(&dir, "u.fq.gz", Encoding::Bgzf)];
        let layout = OutputLayout { unpaired: Some(0), ..OutputLayout::default() };
        run_to(sources(5, 2), 1..=5, &outputs, &config(layout)).unwrap();
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        assert!(std::fs::read(dir.join("u.fq.gz")).unwrap().ends_with(&eof));
    }

    #[test]
    fn read_error_fails_the_run_with_that_error() {
        let dir = scratch_dir("read-error");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let failing = failing_sources(50, 3, 23);
        let err = run_to(failing, 1..=50, &outputs, &config(paired_layout())).unwrap_err();
        assert!(err.to_string().contains("test read failure at spot 23"), "{err:#}");
    }

    #[test]
    fn invalid_spot_fails_the_run() {
        let dir = scratch_dir("invalid");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let mut invalid = sources(10, 2);
        for source in &mut invalid {
            source.spots[4].qualities.pop();
        }
        let err = run_to(invalid, 1..=10, &outputs, &config(paired_layout())).unwrap_err();
        assert!(err.to_string().contains("spot 5"), "{err:#}");
    }

    #[test]
    fn failed_run_writes_no_end_of_file_block() {
        let dir = scratch_dir("failed-eof");
        let outputs = [output(&dir, "u.fq.gz", Encoding::Bgzf)];
        let layout = OutputLayout { unpaired: Some(0), ..OutputLayout::default() };
        assert!(run_to(failing_sources(50, 2, 40), 1..=50, &outputs, &config(layout)).is_err());
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        assert!(!std::fs::read(dir.join("u.fq.gz")).unwrap().ends_with(&eof));
    }

    #[test]
    fn unwritable_output_fails_the_run() {
        let dir = scratch_dir("unwritable");
        let outputs = [Output {
            target: OutputTarget::Path(dir.join("no-such-dir").join("u.fq")),
            encoding: Encoding::Plain,
        }];
        let layout = OutputLayout { unpaired: Some(0), ..OutputLayout::default() };
        let err = run_to(sources(10, 2), 1..=10, &outputs, &config(layout)).unwrap_err();
        assert!(err.to_string().contains("failed to create"), "{err:#}");
    }

    #[test]
    fn unwritable_later_output_fails_the_run_with_its_own_error() {
        let dir = scratch_dir("unwritable-later");
        let outputs = [
            output(&dir, "r1.fq", Encoding::Plain),
            Output {
                target: OutputTarget::Path(dir.join("no-such-dir").join("r2.fq")),
                encoding: Encoding::Plain,
            },
            output(&dir, "u.fq", Encoding::Plain),
        ];
        // Enough batches that dispatch finds the failed writer gone part way through a batch.
        let err =
            run_to(sources(2_000, 2), 1..=2_000, &outputs, &config(paired_layout())).unwrap_err();
        assert!(err.to_string().contains("failed to create"), "{err:#}");
    }

    /// Make a FIFO at `path` with the `mkfifo` utility.
    fn make_fifo(path: &Path) {
        let status = std::process::Command::new("mkfifo").arg(path).status().unwrap();
        assert!(status.success());
    }

    #[test]
    fn fifo_output_is_written_in_place() {
        let dir = scratch_dir("fifo");
        let fifo = dir.join("u.fifo");
        make_fifo(&fifo);
        let reader = {
            let fifo = fifo.clone();
            std::thread::spawn(move || std::fs::read_to_string(fifo).unwrap())
        };
        let layout = OutputLayout { unpaired: Some(0), ..OutputLayout::default() };
        let outputs =
            [Output { target: OutputTarget::Path(fifo.clone()), encoding: Encoding::Plain }];
        run_to(sources(30, 2), 1..=30, &outputs, &config(layout.clone())).unwrap();
        let through_fifo = reader.join().unwrap();

        let file_outputs = [output(&dir, "u.fq", Encoding::Plain)];
        run_to(sources(30, 2), 1..=30, &file_outputs, &config(layout)).unwrap();
        assert_eq!(through_fifo, read_text(&dir.join("u.fq")));
        assert!(std::os::unix::fs::FileTypeExt::is_fifo(
            &std::fs::metadata(&fifo).unwrap().file_type()
        ));
    }

    #[test]
    fn closed_pipe_stops_the_run_early_without_error() {
        let dir = scratch_dir("closed-pipe");
        let fifo = dir.join("u.fifo");
        make_fifo(&fifo);
        let reader = {
            let fifo = fifo.clone();
            std::thread::spawn(move || {
                let mut first_bytes = [0u8; 100];
                File::open(fifo).unwrap().read_exact(&mut first_bytes).unwrap();
            })
        };
        let layout = OutputLayout { unpaired: Some(0), ..OutputLayout::default() };
        let outputs = [Output { target: OutputTarget::Path(fifo), encoding: Encoding::Plain }];
        let summary = run_to(sources(200_000, 2), 1..=200_000, &outputs, &config(layout)).unwrap();
        reader.join().unwrap();
        assert!(summary.stopped_early);
        assert!(summary.counts.spots < 200_000);
    }
}
