//! A BGZF writer that compresses blocks on several threads and writes them in order.

use std::io::{self, Write};
use std::thread::JoinHandle;

use bgzf::{BGZF_BLOCK_SIZE, CompressionLevel, Compressor};
use crossbeam_channel::{Receiver, Sender, bounded};

/// Compression level of BAM and gzip output: libdeflate's level 6, as the default
/// of the zlib-based writers it replaces, so output sizes stay alike.
const COMPRESSION_LEVEL: u8 = 6;

/// Blocks each compressor thread may have queued or in flight, bounding memory.
const BLOCKS_PER_THREAD: usize = 4;

/// An uncompressed block, and where to send it once compressed.
type Job = (Vec<u8>, Sender<io::Result<Vec<u8>>>);

/// A BGZF stream: bytes written are cut into blocks, compressed on `threads` threads, and
/// written to the sink in order by a writer thread.
///
/// [`Self::finish`] ends the stream with the BGZF end-of-file block. A writer dropped
/// without finishing writes no end-of-file block, so its output is detectably truncated.
pub struct ParallelBgzfWriter {
    /// The block being filled.
    buffer: Vec<u8>,
    /// Blocks to compress; `None` once finished or dropped.
    jobs: Option<Sender<Job>>,
    /// Each block's compressed bytes, in the order written; `None` once finished or dropped.
    order: Option<Sender<Receiver<io::Result<Vec<u8>>>>>,
    compressors: Vec<JoinHandle<()>>,
    writer: Option<JoinHandle<io::Result<Box<dyn Write + Send>>>>,
    /// The error that stopped the writer thread, reported by every later call.
    error: Option<io::Error>,
}

impl ParallelBgzfWriter {
    /// Write BGZF to `sink`, compressing on `threads` threads (at least one).
    pub fn new(sink: Box<dyn Write + Send>, threads: usize) -> Self {
        let level = CompressionLevel::new(COMPRESSION_LEVEL).expect("a valid libdeflate level");
        let threads = threads.max(1);
        let (job_tx, job_rx) = bounded::<Job>(threads * BLOCKS_PER_THREAD);
        let (order_tx, order_rx) = bounded(threads * BLOCKS_PER_THREAD);
        let compressors = (0..threads)
            .map(|_| {
                let job_rx = job_rx.clone();
                std::thread::spawn(move || compress_blocks(&job_rx, level))
            })
            .collect();
        let writer = std::thread::spawn(move || write_blocks(&order_rx, sink));
        Self {
            buffer: Vec::with_capacity(BGZF_BLOCK_SIZE),
            jobs: Some(job_tx),
            order: Some(order_tx),
            compressors,
            writer: Some(writer),
            error: None,
        }
    }

    /// Queue the buffered bytes, if any, as one block.
    fn send_block(&mut self) -> io::Result<()> {
        if let Some(error) = &self.error {
            return Err(io::Error::new(error.kind(), error.to_string()));
        }
        if self.buffer.is_empty() {
            return Ok(());
        }
        let block = std::mem::replace(&mut self.buffer, Vec::with_capacity(BGZF_BLOCK_SIZE));
        let (done_tx, done_rx) = bounded(1);
        let sent = match (&self.jobs, &self.order) {
            (Some(jobs), Some(order)) => {
                order.send(done_rx).is_ok() && jobs.send((block, done_tx)).is_ok()
            }
            _ => false,
        };
        if sent {
            return Ok(());
        }
        // The writer thread stops early only on a write error: report that one.
        let error = match self.stop() {
            Err(error) => error,
            Ok(_) => io::Error::other("the BGZF writer stopped"),
        };
        let reported = io::Error::new(error.kind(), error.to_string());
        self.error = Some(error);
        Err(reported)
    }

    /// Write the last block and the end-of-file block, and flush the sink.
    pub fn finish(mut self) -> io::Result<()> {
        self.send_block()?;
        let mut sink = self.stop()?;
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        sink.write_all(&eof)?;
        sink.flush()
    }

    /// Close the channels and wait for every thread, returning the sink, or the writer
    /// thread's error.
    fn stop(&mut self) -> io::Result<Box<dyn Write + Send>> {
        self.jobs = None;
        self.order = None;
        for compressor in self.compressors.drain(..) {
            compressor.join().map_err(|_| io::Error::other("a BGZF compressor panicked"))?;
        }
        match self.writer.take() {
            Some(writer) => {
                writer.join().map_err(|_| io::Error::other("the BGZF writer thread panicked"))?
            }
            None => Err(io::Error::other("BGZF writer already finished")),
        }
    }
}

impl Write for ParallelBgzfWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let taken = buf.len().min(BGZF_BLOCK_SIZE - self.buffer.len());
        self.buffer.extend_from_slice(&buf[..taken]);
        if self.buffer.len() == BGZF_BLOCK_SIZE {
            self.send_block()?;
        }
        Ok(taken)
    }

    /// Queue the buffered bytes as a (short) block. Blocks are written asynchronously, so
    /// this doesn't wait for them to reach the sink; [`Self::finish`] does.
    fn flush(&mut self) -> io::Result<()> {
        self.send_block()
    }
}

impl Drop for ParallelBgzfWriter {
    fn drop(&mut self) {
        if self.writer.is_some() {
            let _ = self.stop();
        }
    }
}

/// Compressor loop: compress each block received, sending it back to its writer slot.
fn compress_blocks(jobs: &Receiver<Job>, level: CompressionLevel) {
    let mut compressor = Compressor::new(level);
    for (block, done) in jobs {
        let mut compressed = Vec::new();
        let result = compressor
            .compress(&block, &mut compressed)
            .map(|()| compressed)
            .map_err(|e| io::Error::other(format!("BGZF compression failed: {e}")));
        // The writer may have gone after a write error; nothing is waiting for this block.
        let _ = done.send(result);
    }
}

/// Writer loop: write each block's compressed bytes in order, returning the sink.
fn write_blocks(
    order: &Receiver<Receiver<io::Result<Vec<u8>>>>,
    mut sink: Box<dyn Write + Send>,
) -> io::Result<Box<dyn Write + Send>> {
    for block in order {
        let compressed = block.recv().map_err(|_| io::Error::other("a BGZF block was lost"))??;
        sink.write_all(&compressed)?;
    }
    Ok(sink)
}

#[cfg(test)]
mod tests {
    use std::io::Read;
    use std::sync::{Arc, Mutex};

    use super::*;

    /// A sink collecting everything written to it.
    #[derive(Clone, Default)]
    struct SharedSink(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedSink {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// A sink that fails every write.
    struct FailingSink;

    impl Write for FailingSink {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("disk full"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn eof_block() -> Vec<u8> {
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        eof
    }

    /// Bytes that span several blocks and don't compress to nothing.
    fn test_bytes(len: usize) -> Vec<u8> {
        (0..len).map(|i| b"ACGTN\n"[(i * 7 + i / 13) % 6]).collect()
    }

    fn decompress(bgzf: &[u8]) -> Vec<u8> {
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(bgzf).read_to_end(&mut plain).unwrap();
        plain
    }

    fn compress_with(data: &[u8], threads: usize, chunk: usize) -> Vec<u8> {
        let sink = SharedSink::default();
        let mut writer = ParallelBgzfWriter::new(Box::new(sink.clone()), threads);
        for piece in data.chunks(chunk) {
            writer.write_all(piece).unwrap();
        }
        writer.finish().unwrap();
        sink.0.lock().unwrap().clone()
    }

    #[test]
    fn output_decompresses_to_the_input_and_ends_with_eof() {
        let data = test_bytes(5 * BGZF_BLOCK_SIZE + 123);
        let compressed = compress_with(&data, 3, 1000);
        assert_eq!(decompress(&compressed), data);
        assert!(compressed.ends_with(&eof_block()));
    }

    #[test]
    fn output_is_identical_at_any_thread_count_or_write_size() {
        let data = test_bytes(3 * BGZF_BLOCK_SIZE + 7);
        let one = compress_with(&data, 1, 1 << 20);
        assert_eq!(compress_with(&data, 8, 1 << 20), one);
        assert_eq!(compress_with(&data, 4, 17), one);
    }

    #[test]
    fn empty_output_is_just_the_eof_block() {
        assert_eq!(compress_with(&[], 2, 1), eof_block());
    }

    #[test]
    fn dropped_writer_writes_no_eof_block() {
        let sink = SharedSink::default();
        {
            let mut writer = ParallelBgzfWriter::new(Box::new(sink.clone()), 2);
            writer.write_all(&test_bytes(2 * BGZF_BLOCK_SIZE)).unwrap();
        }
        let bytes = sink.0.lock().unwrap().clone();
        assert!(!bytes.is_empty());
        assert!(!bytes.ends_with(&eof_block()));
    }

    #[test]
    fn sink_error_is_reported_by_finish() {
        let mut writer = ParallelBgzfWriter::new(Box::new(FailingSink), 2);
        // Writes may or may not see the error, depending on timing; finish always does.
        let _ = writer.write_all(&test_bytes(3 * BGZF_BLOCK_SIZE));
        let err = writer.finish().unwrap_err();
        assert!(err.to_string().contains("disk full"), "{err}");
    }
}
