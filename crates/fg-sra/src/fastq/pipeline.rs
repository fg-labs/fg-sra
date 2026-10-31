//! Converting batches of spots to output bytes on worker threads, and writing them in order.
//!
//! The main thread dispatches fixed-size batches of spot ids. For each batch it hands every
//! output's writer a one-shot receiver, in dispatch order, before handing the batch to the
//! worker pool. A worker reads, routes and formats the batch's spots, compresses each
//! output's bytes into complete BGZF blocks (or leaves them plain), and sends them through
//! the batch's one-shot senders. Each writer waits on its receivers in order, so output is in
//! spot order whichever worker finishes first.
//!
//! How far dispatch runs ahead of the writers is limited by what each output has waiting
//! ([`FlowControl`]), and a batch with nothing for an output stops counting for it once
//! converted, so a reader that takes outputs written in step (such as the two mates) in turns,
//! each turn smaller than the byte limit, can always get to its next records.

use std::fmt;
use std::io::{BufWriter, ErrorKind, Write};
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Condvar, Mutex, MutexGuard, PoisonError};

use anyhow::{Context, Result, anyhow};
use bgzf::{BGZF_BLOCK_SIZE, CompressionLevel, Compressor};
use crossbeam_channel::{Receiver, Sender, bounded, unbounded};

use super::counts::SpotCounts;
use super::format::RecordFormatter;
use super::reader::SpotSource;
use super::spot::{SpotOutcome, SpotRouter};
use super::subsample::Subsampler;
use crate::archive::VDB_THREAD_STACK_BYTES;
use crate::pending_file::PendingFile;
use crate::progress::ProgressLogger;

/// Spots per batch. It is fixed, whatever the thread count, because BGZF blocks end only at
/// batch ends and at every `BGZF_BLOCK_SIZE` bytes within a batch; that keeps output bytes
/// identical at any thread count. A multiple of 8,192, the rows per blob in archives loaded
/// from FASTQ, so that no blob is decoded by two workers.
pub const BATCH_SPOTS: i64 = 16_384;

/// Size of each writer's output buffer.
const WRITE_BUFFER_BYTES: usize = 256 * 1024;

/// Bytes converted for one output but not yet written, below which dispatch never waits. A
/// reader that takes outputs written in step (such as the two mates) in turns is served while
/// each turn, plus a writer's buffer, fits within it: kraken2 takes 10,000 records from one
/// mate's file and then the other's, which fits for reads up to about 3 kb.
pub const MAX_PENDING_BYTES: usize = 64 << 20;

/// Batches dispatched for an output but not yet written, per worker, below which dispatch never
/// waits: enough for every worker to have one batch in hand and one queued. More only holds
/// converted bytes in memory when an output's reader is slower than the workers.
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
    /// Which spots are converted; the rest are skipped unread.
    pub subsampler: Subsampler,
    pub router: SpotRouter,
    pub formatter: RecordFormatter,
    pub layout: OutputLayout,
    /// Whether qualities are read, so spots must carry one per base.
    pub with_qualities: bool,
    /// BGZF compression level, for BGZF outputs.
    pub compression_level: CompressionLevel,
    /// Spots per batch: [`BATCH_SPOTS`], other than in tests.
    pub batch_spots: i64,
    /// Bytes waiting for one output below which dispatch never waits: [`MAX_PENDING_BYTES`],
    /// other than in tests.
    pub max_pending_bytes: usize,
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

/// Whether the run is stopping, and what each output has converted but not yet written, which
/// dispatch waits on.
///
/// Dispatch waits while any one output has both `max_pending_bytes` converted and
/// `min_pending_batches` dispatched but not yet written, so a slow output holds back the run.
/// A batch converted with nothing for an output stops counting for it, so such batches cost
/// that output nothing. The bytes give a reader that takes outputs in turns room to get ahead
/// however sparse the batches; the batches keep every worker busy however large they are.
/// Bytes are known only once a batch is converted, so the batches being converted can take an
/// output past the byte limit.
///
/// No limit serves a reader that reads one output only after finishing another, or that takes
/// outputs not written in step (such as pairs and unpaired reads) in turns: with enough spots,
/// an output it isn't reading reaches the limit first.
struct FlowControl {
    stop: AtomicBool,
    max_pending_bytes: usize,
    min_pending_batches: usize,
    pending: Mutex<Pending>,
    /// Signalled when a count dispatch waits on falls, a writer exits, or the run stops.
    changed: Condvar,
}

/// What dispatch waits on.
struct Pending {
    /// Bytes converted but not yet written, per output.
    bytes: Vec<usize>,
    /// Batches dispatched but not yet written, per output, less those converted with nothing
    /// for it.
    batches: Vec<usize>,
    /// A writer has exited. Before dispatch finishes that is always a failure, which dispatch
    /// finds on its next send to the writer, so it waits no longer.
    writer_exited: bool,
}

impl FlowControl {
    /// Flow control for `num_outputs` outputs, waiting once an output has `max_pending_bytes`
    /// and `min_pending_batches` (each at least 1) waiting.
    fn new(num_outputs: usize, max_pending_bytes: usize, min_pending_batches: usize) -> Self {
        Self {
            stop: AtomicBool::new(false),
            max_pending_bytes: max_pending_bytes.max(1),
            min_pending_batches: min_pending_batches.max(1),
            pending: Mutex::new(Pending {
                bytes: vec![0; num_outputs],
                batches: vec![0; num_outputs],
                writer_exited: false,
            }),
            changed: Condvar::new(),
        }
    }

    fn stopping(&self) -> bool {
        self.stop.load(Ordering::Relaxed)
    }

    /// Stop the run, waking dispatch if it is waiting.
    fn stop(&self) {
        // Set under the lock, so dispatch can't check the flag, miss this, and then wait.
        let pending = self.lock();
        self.stop.store(true, Ordering::Relaxed);
        drop(pending);
        self.changed.notify_all();
    }

    /// Count a converted batch's bytes for each output, before they are sent to the writers.
    /// Count a batch against every output, before it is handed to the workers.
    fn count_dispatched(&self) {
        for batches in &mut self.lock().batches {
            *batches += 1;
        }
    }

    /// Count a converted batch's bytes for each output, before they are sent to the writers. A
    /// batch with nothing for an output stops counting for it.
    fn count_converted(&self, batch: &[Vec<u8>]) {
        let mut pending = self.lock();
        let Pending { bytes, batches, .. } = &mut *pending;
        for ((output_bytes, output_batches), batch_bytes) in
            bytes.iter_mut().zip(batches.iter_mut()).zip(batch)
        {
            if batch_bytes.is_empty() {
                debug_assert!(*output_batches > 0, "more batches converted than dispatched");
                *output_batches = output_batches.saturating_sub(1);
            } else {
                *output_bytes += batch_bytes.len();
            }
        }
        drop(pending);
        self.changed.notify_all();
    }

    /// Count a batch of `len` bytes of `output` as written. One with no bytes stopped counting
    /// when it was converted.
    fn count_written(&self, output: usize, len: usize) {
        if len == 0 {
            return;
        }
        let mut pending = self.lock();
        debug_assert!(pending.bytes[output] >= len, "more bytes written than converted");
        debug_assert!(pending.batches[output] > 0, "more batches written than dispatched");
        pending.bytes[output] = pending.bytes[output].saturating_sub(len);
        pending.batches[output] = pending.batches[output].saturating_sub(1);
        drop(pending);
        self.changed.notify_all();
    }

    /// Note that a writer has exited, for whatever reason.
    fn writer_exited(&self) {
        self.lock().writer_exited = true;
        self.changed.notify_all();
    }

    /// Wait until no output has both `max_pending_bytes` and `min_pending_batches` waiting, a
    /// writer has exited, or the run is stopping.
    fn wait_for_room(&self) {
        let pending = self.lock();
        let _pending = self
            .changed
            .wait_while(pending, |pending| {
                !self.stopping()
                    && !pending.writer_exited
                    && pending.bytes.iter().zip(&pending.batches).any(|(&bytes, &batches)| {
                        bytes >= self.max_pending_bytes && batches >= self.min_pending_batches
                    })
            })
            .unwrap_or_else(PoisonError::into_inner);
    }

    /// What dispatch waits on. A poisoned lock is used as it is: each value changes in one
    /// step, so a panic can't leave one half-changed.
    fn lock(&self) -> MutexGuard<'_, Pending> {
        self.pending.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

/// Notes that a writer has exited when dropped, however the writer exits.
struct WriterExit<'a>(&'a FlowControl);

impl Drop for WriterExit<'_> {
    fn drop(&mut self) {
        self.0.writer_exited();
    }
}

/// Convert the spots in `spots` to `outputs`, with one worker thread per source.
///
/// On failure, the first real error is returned: a worker's (reading, validating or
/// compressing) in preference to a writer's. If an output's reader goes away the run stops
/// early without error, and the summary says so. BGZF outputs get their end-of-file block
/// only when every spot was converted.
///
/// Also returns the outputs written under a temporary name, to rename into place with
/// [`PendingFile::commit_all`] once the run is known to be good; dropping them removes them.
pub fn run<S: SpotSource + Send>(
    sources: Vec<S>,
    spots: RangeInclusive<i64>,
    outputs: &[Output],
    config: &PipelineConfig,
    progress: &ProgressLogger,
) -> Result<(PipelineSummary, Vec<PendingFile>)> {
    let threads = sources.len().max(1);
    let flow =
        FlowControl::new(outputs.len(), config.max_pending_bytes, threads * BATCHES_PER_WORKER);
    let pipe_closed = AtomicBool::new(false);

    std::thread::scope(|scope| {
        let (job_tx, job_rx) = bounded::<Job>(threads);
        let mut order_txs = Vec::with_capacity(outputs.len());
        let mut writers = Vec::with_capacity(outputs.len());
        for (index, output) in outputs.iter().enumerate() {
            // Unbounded: `flow` limits how far dispatch runs ahead, by bytes.
            let (order_tx, order_rx) = unbounded();
            order_txs.push(order_tx);
            let (flow, pipe_closed) = (&flow, &pipe_closed);
            writers.push(
                scope.spawn(move || write_output(&order_rx, index, output, flow, pipe_closed)),
            );
        }
        let workers = sources
            .into_iter()
            .map(|source| {
                let (job_rx, flow) = (job_rx.clone(), &flow);
                std::thread::Builder::new()
                    .stack_size(VDB_THREAD_STACK_BYTES)
                    .spawn_scoped(scope, move || {
                        work(source, &job_rx, outputs, config, flow, progress)
                    })
            })
            .collect::<std::io::Result<Vec<_>>>()
            .context("failed to start a worker thread")?;
        drop(job_rx);

        dispatch(spots, config.batch_spots, &job_tx, &order_txs, &flow);
        drop(job_tx);
        drop(order_txs);

        let mut summary = PipelineSummary::default();
        let mut pending = Vec::new();
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
            match writer.join().unwrap_or_else(|_| Err(anyhow!("a writer thread panicked"))) {
                Ok(file) => pending.extend(file),
                Err(e) => keep_first_real_error(&mut first_error, e),
            }
        }
        if let Some(e) = first_error {
            return Err(e);
        }
        summary.stopped_early = pipe_closed.load(Ordering::Relaxed);
        Ok((summary, pending))
    })
}

/// Keep `error` as the run's error unless one is kept already or it is only [`Abandoned`].
fn keep_first_real_error(first: &mut Option<anyhow::Error>, error: anyhow::Error) {
    if first.is_none() && !error.is::<Abandoned>() {
        *first = Some(error);
    }
}

/// Hand out `spots` in batches until done or told to stop, waiting while any output has too
/// many bytes waiting. Each writer is given its receiver before the batch goes to a worker, so
/// no batch's bytes can arrive unexpected.
fn dispatch(
    spots: RangeInclusive<i64>,
    batch_spots: i64,
    job_tx: &Sender<Job>,
    order_txs: &[Sender<Receiver<Result<Vec<u8>>>>],
    flow: &FlowControl,
) {
    let (mut first, end) = spots.into_inner();
    while first <= end {
        flow.wait_for_room();
        if flow.stopping() {
            return;
        }
        let last = end.min(first + batch_spots.max(1) - 1);
        let mut results = Vec::with_capacity(order_txs.len());
        for order_tx in order_txs {
            let (result_tx, result_rx) = bounded(1);
            if order_tx.send(result_rx).is_err() {
                // The writer has failed; its error is reported when it is joined. The writers
                // already given this batch skip it rather than report it as never converted.
                flow.stop();
                Job { first, last, results }.abandon();
                return;
            }
            results.push(result_tx);
        }
        flow.count_dispatched();
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
    flow: &FlowControl,
    progress: &ProgressLogger,
) -> Result<WorkerTally> {
    let num_outputs = outputs.len();
    let mut compressors: Vec<_> =
        (0..num_outputs).map(|_| Compressor::new(config.compression_level)).collect();
    let mut buffers = vec![Vec::new(); num_outputs];
    let mut tally = WorkerTally::default();
    let mut failure = None;

    for job in jobs {
        if failure.is_some() || flow.stopping() {
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
                flow.count_converted(&encoded);
                for (result, bytes) in job.results.into_iter().zip(encoded) {
                    // A writer that has already exited no longer needs the bytes.
                    let _ = result.send(Ok(bytes));
                }
                progress.record((job.last - job.first + 1) as u64);
            }
            Err(e) => {
                flow.stop();
                job.abandon();
                failure = Some(e);
            }
        }
    }
    failure.map_or(Ok(tally), Err)
}

/// Read, route and format the spots of `ids` that the subsampler keeps into one buffer per
/// output.
fn convert_batch<S: SpotSource>(
    source: &mut S,
    ids: RangeInclusive<i64>,
    config: &PipelineConfig,
    buffers: &mut [Vec<u8>],
) -> Result<WorkerTally> {
    let layout = &config.layout;
    let formatter = &config.formatter;
    let mut tally = WorkerTally::default();
    let kept: Vec<i64> = ids.clone().filter(|&id| config.subsampler.keeps(id)).collect();
    let batch_len = (ids.end() - ids.start() + 1) as u64;
    tally.counts.spots_skipped_by_subsampling = batch_len - kept.len() as u64;
    source.start_batch(&kept)?;
    for &id in &kept {
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

/// Writer loop for one output: write each batch's bytes in dispatch order, returning the
/// output's pending file (if written under a temporary name) for the caller to commit.
///
/// A closed pipe (the output's reader went away) stops the whole run: this writer notes it,
/// then drains the batches already dispatched without writing them. The BGZF end-of-file
/// block is written only if the run wasn't stopped, so an incomplete output is detectably
/// truncated.
fn write_output(
    order: &Receiver<Receiver<Result<Vec<u8>>>>,
    index: usize,
    output: &Output,
    flow: &FlowControl,
    pipe_closed: &AtomicBool,
) -> Result<Option<PendingFile>> {
    let _exit = WriterExit(flow);
    let name = output.target.describe();
    let (sink, pending): (Box<dyn Write>, _) = match &output.target {
        OutputTarget::Stdout => (Box::new(std::io::stdout().lock()), None),
        OutputTarget::Path(path) => {
            let (file, pending) = PendingFile::create(path)?;
            (Box::new(file), pending)
        }
    };
    let mut writer = BufWriter::with_capacity(WRITE_BUFFER_BYTES, sink);
    let mut closed = false;
    let note_closed = |closed: &mut bool| {
        *closed = true;
        pipe_closed.store(true, Ordering::Relaxed);
        flow.stop();
    };

    for batch in order {
        let bytes = match batch.recv() {
            Ok(Ok(bytes)) => bytes,
            Ok(Err(e)) if e.is::<Abandoned>() => continue,
            Ok(Err(e)) => return Err(e),
            Err(_) => return Err(anyhow!("a worker exited without converting a batch")),
        };
        if !closed && let Err(e) = writer.write_all(&bytes) {
            if e.kind() != ErrorKind::BrokenPipe {
                return Err(e).with_context(|| format!("failed to write {name}"));
            }
            note_closed(&mut closed);
        }
        flow.count_written(index, bytes.len());
    }
    if closed || flow.stopping() {
        return Ok(pending);
    }
    if output.encoding == Encoding::Bgzf {
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        if let Err(e) = writer.write_all(&eof) {
            if e.kind() != ErrorKind::BrokenPipe {
                return Err(e).with_context(|| format!("failed to write {name}"));
            }
            note_closed(&mut closed);
            return Ok(pending);
        }
    }
    match writer.flush() {
        Err(e) if e.kind() == ErrorKind::BrokenPipe => {
            note_closed(&mut closed);
            Ok(pending)
        }
        result => result.map(|()| pending).with_context(|| format!("failed to write {name}")),
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::Read;
    use std::sync::atomic::AtomicI64;
    use std::sync::{Arc, Mutex};

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

    /// A source whose every read first uses 4 MiB of stack, twice Rust's default for a thread;
    /// libncbi-vdb can need more. (Debug builds may copy the array, so it is kept well under
    /// the 16 MiB given.)
    struct DeepStackSource(TestSource);

    impl SpotSource for DeepStackSource {
        // A large stack frame is the point.
        #[allow(clippy::large_stack_arrays)]
        fn read(&mut self, id: i64) -> Result<Spot<'_>> {
            let stack = std::hint::black_box([1u8; 4 << 20]);
            assert_eq!(stack[stack.len() - 1], 1);
            self.0.read(id)
        }
    }

    /// A source that logs each batch it is told of, and fails a read outside that batch.
    struct BatchLoggingSource {
        source: TestSource,
        batch: Option<Vec<i64>>,
        log: Arc<Mutex<Vec<Vec<i64>>>>,
    }

    impl SpotSource for BatchLoggingSource {
        fn start_batch(&mut self, ids: &[i64]) -> Result<()> {
            self.log.lock().unwrap().push(ids.to_vec());
            self.batch = Some(ids.to_vec());
            Ok(())
        }

        fn read(&mut self, id: i64) -> Result<Spot<'_>> {
            if !self.batch.as_ref().is_some_and(|batch| batch.contains(&id)) {
                anyhow::bail!("spot {id} read outside the batch started, {:?}", self.batch);
            }
            self.source.read(id)
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
            subsampler: Subsampler::new(1.0, 42),
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
            // Small, so that dispatch waits on the writers often.
            max_pending_bytes: 1024,
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
        let (summary, pending) = run(sources, spots, outputs, config, &ProgressLogger::new(0, 0))?;
        PendingFile::commit_all(pending)?;
        Ok(summary)
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
    fn each_batch_is_started_before_its_spots_are_read() {
        let dir = scratch_dir("batches");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let log = Arc::default();
        let sources: Vec<_> = sources(10, 3)
            .into_iter()
            .map(|source| BatchLoggingSource { source, batch: None, log: Arc::clone(&log) })
            .collect();
        run_to(sources, 1..=10, &outputs, &config(paired_layout())).unwrap();
        let mut batches = log.lock().unwrap().clone();
        batches.sort();
        assert_eq!(batches, [vec![1, 2, 3, 4], vec![5, 6, 7, 8], vec![9, 10]]);
    }

    /// `config` for `layout`, keeping `fraction` of spots by seed 42.
    fn subsampled_config(layout: OutputLayout, fraction: f64) -> PipelineConfig {
        PipelineConfig { subsampler: Subsampler::new(fraction, 42), ..config(layout) }
    }

    /// The spot ids of `fastq`'s records, named `$ac.$si/$ri`.
    fn spot_ids(fastq: &str) -> Vec<i64> {
        names(fastq)
            .iter()
            .map(|name| name.trim_start_matches("@T.").split('/').next().unwrap().parse().unwrap())
            .collect()
    }

    #[test]
    fn subsampling_keeps_or_skips_each_spot_in_every_output_together() {
        let dir = scratch_dir("subsample-mates");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let config = subsampled_config(paired_layout(), 0.5);
        run_to(sources(90, 2), 1..=90, &outputs, &config).unwrap();

        let first_mates = spot_ids(&read_text(&dir.join("r1.fq")));
        assert_eq!(first_mates, spot_ids(&read_text(&dir.join("r2.fq"))));
        let mut written = [first_mates, spot_ids(&read_text(&dir.join("u.fq")))].concat();
        written.sort_unstable();
        let kept: Vec<i64> = (1..=90).filter(|&id| config.subsampler.keeps(id)).collect();
        assert_eq!(written, kept);
        assert!(!kept.is_empty() && kept.len() < 90, "{kept:?}");
    }

    #[test]
    fn spots_skipped_by_subsampling_are_never_read_but_are_counted() {
        let dir = scratch_dir("subsample-unread");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let log = Arc::default();
        let sources: Vec<_> = sources(40, 3)
            .into_iter()
            .map(|source| BatchLoggingSource { source, batch: None, log: Arc::clone(&log) })
            .collect();
        let config = subsampled_config(paired_layout(), 0.3);
        let summary = run_to(sources, 1..=40, &outputs, &config).unwrap();

        let mut started: Vec<i64> = log.lock().unwrap().concat();
        started.sort_unstable();
        let kept: Vec<i64> = (1..=40).filter(|&id| config.subsampler.keeps(id)).collect();
        assert_eq!(started, kept);
        assert_eq!(summary.counts.spots, kept.len() as u64);
        assert_eq!(summary.counts.spots_skipped_by_subsampling, 40 - kept.len() as u64);
    }

    #[test]
    fn subsampled_output_is_identical_at_one_and_four_threads() {
        let convert = |threads| {
            let dir = scratch_dir(&format!("subsample-threads-{threads}"));
            let outputs =
                ["r1.fq.gz", "r2.fq.gz", "u.fq.gz"].map(|n| output(&dir, n, Encoding::Bgzf));
            let config = subsampled_config(paired_layout(), 0.2);
            run_to(sources(500, threads), 1..=500, &outputs, &config).unwrap();
            outputs.map(|output| match output.target {
                OutputTarget::Path(path) => std::fs::read(path).unwrap(),
                OutputTarget::Stdout => unreachable!(),
            })
        };
        assert!(convert(1) == convert(4), "subsampled output differs between 1 and 4 threads");
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
    fn failed_run_leaves_existing_outputs_untouched() {
        let dir = scratch_dir("failed-existing");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        for name in ["r1.fq", "u.fq"] {
            std::fs::write(dir.join(name), b"old\n").unwrap();
        }
        assert!(
            run_to(failing_sources(50, 2, 40), 1..=50, &outputs, &config(paired_layout())).is_err()
        );
        assert_eq!(read_text(&dir.join("r1.fq")), "old\n");
        assert_eq!(read_text(&dir.join("u.fq")), "old\n");
        let mut left: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        left.sort();
        assert_eq!(left, ["r1.fq", "u.fq"], "no new or temporary files are left behind");
    }

    /// An output written in place (here through a symbolic link) can't be held back, so a
    /// failed run leaves it without the BGZF end-of-file block, detectably truncated.
    #[cfg(unix)]
    #[test]
    fn failed_run_writes_no_end_of_file_block() {
        let dir = scratch_dir("failed-eof");
        let target = dir.join("target.fq.gz");
        std::os::unix::fs::symlink(&target, dir.join("u.fq.gz")).unwrap();
        let outputs = [output(&dir, "u.fq.gz", Encoding::Bgzf)];
        let layout = OutputLayout { unpaired: Some(0), ..OutputLayout::default() };
        assert!(run_to(failing_sources(50, 2, 40), 1..=50, &outputs, &config(layout)).is_err());
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        assert!(!std::fs::read(&target).unwrap().ends_with(&eof));
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

    /// A spot of biological reads of any length, each base of quality 30 ([`TestSpot::new`]
    /// takes only short reads).
    fn long_read_spot(id: i64, reads: &[&str]) -> TestSpot {
        let mut read_starts = Vec::new();
        let mut bases = Vec::new();
        for read in reads {
            read_starts.push(bases.len() as i32);
            bases.extend_from_slice(read.as_bytes());
        }
        TestSpot {
            id,
            name: Vec::new(),
            spot_group: Vec::new(),
            qualities: vec![30; bases.len()],
            bases,
            read_starts,
            read_lens: reads.iter().map(|read| read.len() as u32).collect(),
            read_types: vec![B; reads.len()],
            read_filters: vec![0; reads.len()],
        }
    }

    /// The names of the next `count` records of `fastq`, or of those left if fewer.
    fn next_names(fastq: impl std::io::BufRead, count: usize) -> Vec<String> {
        let mut names = Vec::with_capacity(count);
        let mut lines = fastq.lines();
        while names.len() < count {
            let Some(name) = lines.next() else { break };
            names.push(name.unwrap());
            for _ in 0..3 {
                lines.next().unwrap().unwrap();
            }
        }
        names
    }

    #[test]
    fn reader_taking_mates_in_turns_gets_past_a_long_run_of_single_read_spots() {
        // Records of about 1 KB, so each of the reader's turns (as kraken2 -P takes them), about
        // 2 MB, is far more than a writer's buffer and its FIFO's hold, but well within the
        // byte limit; and several batches of pairs hold more than a turn.
        const READ_LEN: usize = 500;
        const RECORDS_PER_TURN: usize = 2_000;
        const SPOTS_PER_BATCH: i64 = 1_024;
        const PENDING_LIMIT: usize = 4 << 20;
        // The reader's second turn on r1 straddles the start of 8 batches of single-read
        // spots, which give r1 nothing until the pairs after them.
        const PAIRED_BEFORE: i64 = 3_000;
        const SINGLE: i64 = 8 * SPOTS_PER_BATCH;
        const PAIRED_AFTER: i64 = 3_000;
        let total = PAIRED_BEFORE + SINGLE + PAIRED_AFTER;
        let is_paired = |id: i64| id <= PAIRED_BEFORE || id > PAIRED_BEFORE + SINGLE;

        let dir = scratch_dir("mates-in-turns");
        let (r1, r2) = (dir.join("r1.fifo"), dir.join("r2.fifo"));
        make_fifo(&r1);
        make_fifo(&r2);
        let outputs = [&r1, &r2].map(|path| Output {
            target: OutputTarget::Path(path.clone()),
            encoding: Encoding::Plain,
        });
        let layout =
            OutputLayout { first_mates: Some(0), second_mates: Some(1), ..OutputLayout::default() };
        let read = "ACGT".repeat(READ_LEN / 4);
        let sources: Vec<_> = (0..2)
            .map(|_| TestSource {
                spots: (1..=total)
                    .map(|id| {
                        let reads: &[&str] = if is_paired(id) { &[&read, &read] } else { &[&read] };
                        long_read_spot(id, reads)
                    })
                    .collect(),
                fail_at: None,
            })
            .collect();

        let fifos = [r1.clone(), r2.clone()];
        let reader = std::thread::spawn(move || {
            let mut r1 = std::io::BufReader::new(File::open(r1).unwrap());
            let mut r2 = std::io::BufReader::new(File::open(r2).unwrap());
            let (mut first_mates, mut second_mates) = (Vec::new(), Vec::new());
            loop {
                let turn = next_names(&mut r1, RECORDS_PER_TURN);
                second_mates.extend(next_names(&mut r2, turn.len()));
                let last_turn = turn.len() < RECORDS_PER_TURN;
                first_mates.extend(turn);
                if last_turn {
                    return (first_mates, second_mates);
                }
            }
        });
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        let config = PipelineConfig {
            batch_spots: SPOTS_PER_BATCH,
            max_pending_bytes: PENDING_LIMIT,
            ..config(layout)
        };
        std::thread::spawn(move || {
            let _ = done_tx.send(run_to(sources, 1..=total, &outputs, &config));
        });
        let Ok(result) = done_rx.recv_timeout(std::time::Duration::from_secs(20)) else {
            // Drain both FIFOs so the pipeline's threads can finish.
            for fifo in fifos {
                std::thread::spawn(move || {
                    std::io::copy(&mut File::open(fifo).unwrap(), &mut std::io::sink())
                });
            }
            panic!("the conversion deadlocked with a reader taking mates in turns");
        };
        let summary = result.unwrap();

        let (first_mates, second_mates) = reader.join().unwrap();
        let expected: Vec<i64> = (1..=total).filter(|&id| is_paired(id)).collect();
        let spot_ids = |names: &[String], suffix: &str| -> Vec<i64> {
            names
                .iter()
                .map(|name| name.strip_prefix("@T.").unwrap().strip_suffix(suffix).unwrap())
                .map(|id| id.parse().unwrap())
                .collect()
        };
        assert_eq!(spot_ids(&first_mates, "/1"), expected);
        assert_eq!(spot_ids(&second_mates, "/2"), expected);
        assert_eq!(summary.counts.pairs_written, expected.len() as u64);
    }

    /// A source that panics reading spot `panic_at`.
    struct PanickingSource {
        source: TestSource,
        panic_at: i64,
    }

    impl SpotSource for PanickingSource {
        fn read(&mut self, id: i64) -> Result<Spot<'_>> {
            assert!(id != self.panic_at, "test panic at spot {id}");
            self.source.read(id)
        }
    }

    #[test]
    fn worker_panic_fails_the_run_rather_than_hanging() {
        let dir = scratch_dir("worker-panic");
        let outputs = ["r1.fq", "r2.fq", "u.fq"].map(|n| output(&dir, n, Encoding::Plain));
        let sources: Vec<_> = sources(2_000, 4)
            .into_iter()
            .map(|source| PanickingSource { source, panic_at: 203 })
            .collect();
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let result = run_to(sources, 1..=2_000, &outputs, &config(paired_layout()));
            let _ = done_tx.send(result.map(|_| ()));
        });
        let err = done_rx
            .recv_timeout(std::time::Duration::from_secs(20))
            .expect("the run hung after a worker panicked")
            .unwrap_err();
        assert!(err.to_string().contains("panicked"), "{err:#}");
    }

    /// A source that notes the highest spot id read.
    struct HighestIdSource {
        source: TestSource,
        highest: Arc<AtomicI64>,
    }

    impl SpotSource for HighestIdSource {
        fn read(&mut self, id: i64) -> Result<Spot<'_>> {
            self.highest.fetch_max(id, Ordering::Relaxed);
            self.source.read(id)
        }
    }

    /// Order queues are unbounded, so only the byte limit stops conversion running ahead of an
    /// output's reader.
    #[test]
    fn output_not_being_read_holds_back_conversion() {
        const SPOTS: i64 = 5_000;
        let dir = scratch_dir("held-back");
        let fifo = dir.join("u.fifo");
        make_fifo(&fifo);
        let outputs =
            [Output { target: OutputTarget::Path(fifo.clone()), encoding: Encoding::Plain }];
        let layout = OutputLayout { unpaired: Some(0), ..OutputLayout::default() };
        let highest = Arc::new(AtomicI64::new(0));
        let read = "ACGT".repeat(125);
        let sources: Vec<_> = (0..2)
            .map(|_| HighestIdSource {
                source: TestSource {
                    spots: (1..=SPOTS).map(|id| long_read_spot(id, &[&read])).collect(),
                    fail_at: None,
                },
                highest: Arc::clone(&highest),
            })
            .collect();

        let reader = {
            let highest = Arc::clone(&highest);
            std::thread::spawn(move || {
                let mut fastq = File::open(fifo).unwrap();
                // Long enough for every spot to be converted if nothing held conversion back.
                std::thread::sleep(std::time::Duration::from_millis(500));
                let highest_before_reading = highest.load(Ordering::Relaxed);
                let mut text = String::new();
                fastq.read_to_string(&mut text).unwrap();
                (highest_before_reading, text)
            })
        };
        run_to(sources, 1..=SPOTS, &outputs, &config(layout)).unwrap();
        let (highest_before_reading, text) = reader.join().unwrap();

        // The writer's 256 KiB buffer and the FIFO's 16-64 KiB hold 270-320 records of 1 KB;
        // a few more are waiting or being converted.
        assert!(highest_before_reading < 1_000, "{highest_before_reading} spots were read");
        assert_eq!(names(&text).len(), SPOTS as usize);
    }
}
