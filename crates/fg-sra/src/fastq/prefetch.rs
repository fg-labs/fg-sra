//! Prefetching the aligned reads whose alignments lie far from the rest of their batch's.
//!
//! Spots are converted in spot order, but an aligned archive stores its alignments in
//! reference order. Most of a batch's alignments lie together in the alignment table, so each
//! blob there is decoded once for many reads. A few percent do not (mates aligned far apart,
//! spots numbered out of alignment order), and each of those decodes a blob that no other read
//! of its batch needs. They are found before conversion from `PRIMARY_ALIGNMENT_ID` alone and
//! rebuilt in alignment order, where neighbouring ones share their blobs.

use std::ops::RangeInclusive;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use anyhow::{Context, Result, anyhow};
use fg_sra_vdb::database::VTable;

use super::reader::{AlignmentReader, References};
use crate::archive::VDB_THREAD_STACK_BYTES;

/// Alignment rows per block when judging whether an alignment is near the rest of its batch's:
/// about one blob of the alignment table's columns.
const BLOCK_ROWS_SHIFT: u32 = 10;

/// A batch's alignments in a block holding fewer of its alignments than this are prefetched.
const MIN_BLOCK_ALIGNMENTS: usize = 8;

/// Alignments rebuilt per unit of prefetch work.
const CHUNK_ALIGNMENTS: usize = 4096;

/// At most this many alignments are prefetched, bounding the memory they take (about 2.4 GB at
/// 150 bases each). Archives whose spots are far out of alignment order would otherwise hold
/// most of their reads in memory.
const MAX_PREFETCHED: usize = 16_000_000;

/// Aligned reads rebuilt ahead of conversion, looked up by alignment id.
#[derive(Debug, Default)]
pub struct PrefetchedReads {
    /// Bit `id` is set for each prefetched alignment id, so most lookups need no search.
    present: Vec<u64>,
    /// Prefetched alignment ids, ascending.
    ids: Vec<i64>,
    /// `bases[starts[i]..starts[i + 1]]` is alignment `ids[i]`'s read, in reference orientation.
    starts: Vec<usize>,
    bases: Vec<u8>,
}

impl PrefetchedReads {
    /// Find the alignments of `spots` that lie far from the rest of their batch's, taking batches
    /// of `batch_spots` from the range's start as the conversion does, and rebuild their reads from
    /// `alignments` and `references`, on `threads` threads.
    ///
    /// `table` is the reads table, read for `PRIMARY_ALIGNMENT_ID` only.
    pub fn far_alignments(
        table: &VTable,
        alignments: &VTable,
        references: References<'_>,
        spots: RangeInclusive<i64>,
        batch_spots: i64,
        threads: usize,
    ) -> Result<Self> {
        let mut far = find_far_alignments(table, spots, batch_spots, threads)?;
        far.truncate(MAX_PREFETCHED);
        let readers = (0..threads.max(1))
            .map(|_| AlignmentReader::new(alignments, references))
            .collect::<Result<Vec<_>>>()?;
        let chunks = parallel_chunks(
            &far.chunks(CHUNK_ALIGNMENTS).collect::<Vec<_>>(),
            readers,
            |reader, ids| {
                let mut chunk = Chunk {
                    ids: ids.to_vec(),
                    starts: Vec::with_capacity(ids.len() + 1),
                    bases: Vec::new(),
                };
                let mut read = Vec::new();
                for &id in ids {
                    reader.restore(id, &mut read)?;
                    chunk.starts.push(chunk.bases.len());
                    chunk.bases.extend_from_slice(&read);
                }
                chunk.starts.push(chunk.bases.len());
                Ok(chunk)
            },
        )?;
        Ok(Self::from_chunks(chunks))
    }

    /// Reads from chunks in ascending id order, each chunk's ids ascending.
    fn from_chunks(chunks: Vec<Chunk>) -> Self {
        let mut reads = Self::default();
        let count = chunks.iter().map(|c| c.ids.len()).sum();
        reads.ids.reserve(count);
        reads.starts.reserve(count + 1);
        reads.bases.reserve(chunks.iter().map(|c| c.bases.len()).sum());
        for chunk in chunks {
            let offset = reads.bases.len();
            reads.ids.extend_from_slice(&chunk.ids);
            reads.starts.extend(chunk.starts[..chunk.ids.len()].iter().map(|s| s + offset));
            reads.bases.extend_from_slice(&chunk.bases);
        }
        reads.starts.push(reads.bases.len());
        let max_id = reads.ids.last().map_or(0, |&id| usize::try_from(id).unwrap_or(0));
        reads.present = vec![0; max_id / 64 + 1];
        for &id in &reads.ids {
            let bit = usize::try_from(id).unwrap_or(0);
            reads.present[bit / 64] |= 1 << (bit % 64);
        }
        reads
    }

    /// Alignment `id`'s read, in reference orientation, if it was prefetched.
    #[must_use]
    pub fn get(&self, id: i64) -> Option<&[u8]> {
        let bit = usize::try_from(id).ok()?;
        if self.present.get(bit / 64)? & (1 << (bit % 64)) == 0 {
            return None;
        }
        let index = self.ids.binary_search(&id).ok()?;
        Some(&self.bases[self.starts[index]..self.starts[index + 1]])
    }

    /// How many reads were prefetched.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Bases across every prefetched read.
    #[must_use]
    pub fn total_bases(&self) -> usize {
        self.bases.len()
    }
}

/// Consecutive prefetched reads: `bases[starts[i]..starts[i + 1]]` is alignment `ids[i]`'s.
struct Chunk {
    ids: Vec<i64>,
    starts: Vec<usize>,
    bases: Vec<u8>,
}

/// The alignments of every batch of `spots` that [`far_in_batch`] picks, ascending.
fn find_far_alignments(
    table: &VTable,
    spots: RangeInclusive<i64>,
    batch_spots: i64,
    threads: usize,
) -> Result<Vec<i64>> {
    let (first, last) = spots.into_inner();
    let batch_spots = batch_spots.max(1);
    let step = usize::try_from(batch_spots).unwrap_or(1);
    let batches: Vec<RangeInclusive<i64>> = (first..=last)
        .step_by(step)
        .map(|start| start..=last.min(start + batch_spots - 1))
        .collect();
    let open_cursor = || -> Result<_> {
        let cursor = table.create_cursor_read()?;
        let column = cursor
            .add_column("(I64)PRIMARY_ALIGNMENT_ID")
            .context("failed to read column PRIMARY_ALIGNMENT_ID")?;
        cursor.open().context("failed to open a cursor on the reads table")?;
        Ok((cursor, column, Vec::new()))
    };
    let cursors = (0..threads.max(1)).map(|_| open_cursor()).collect::<Result<Vec<_>>>()?;
    let per_batch =
        parallel_chunks(&batches.iter().collect::<Vec<_>>(), cursors, |state, batch| {
            let (cursor, column, cell) = state;
            let mut ids = Vec::new();
            for spot in batch.clone() {
                cursor
                    .read_i64_slice_into(spot, *column, cell)
                    .with_context(|| format!("failed to read spot {spot}'s alignment ids"))?;
                ids.extend_from_slice(cell);
            }
            Ok(far_in_batch(&mut ids))
        })?;
    let mut far: Vec<i64> = per_batch.into_iter().flatten().collect();
    far.sort_unstable();
    far.dedup();
    Ok(far)
}

/// The alignments among one batch's `ids` (any order, 0 for none; sorted in place) that share
/// their block of alignment rows with fewer than `MIN_BLOCK_ALIGNMENTS` of the batch's
/// alignments, ascending.
fn far_in_batch(ids: &mut Vec<i64>) -> Vec<i64> {
    ids.retain(|&id| id > 0);
    ids.sort_unstable();
    ids.chunk_by(|a, b| a >> BLOCK_ROWS_SHIFT == b >> BLOCK_ROWS_SHIFT)
        .filter(|block| block.len() < MIN_BLOCK_ALIGNMENTS)
        .flatten()
        .copied()
        .collect()
}

/// Run `work` over `items` on one thread per state in `states`, each thread working with its
/// own state, and return the results in item order. The first error stops the rest.
fn parallel_chunks<T: Copy + Sync, S: Send, R: Send>(
    items: &[T],
    states: Vec<S>,
    work: impl Fn(&mut S, T) -> Result<R> + Sync,
) -> Result<Vec<R>> {
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let results: Mutex<Vec<Option<R>>> = Mutex::new((0..items.len()).map(|_| None).collect());
    std::thread::scope(|scope| -> Result<()> {
        let handles = states
            .into_iter()
            .map(|mut state| {
                let (next, failed, results, work) = (&next, &failed, &results, &work);
                std::thread::Builder::new().stack_size(VDB_THREAD_STACK_BYTES).spawn_scoped(
                    scope,
                    move || -> Result<()> {
                        loop {
                            let index = next.fetch_add(1, Ordering::Relaxed);
                            if index >= items.len() || failed.load(Ordering::Relaxed) {
                                return Ok(());
                            }
                            match work(&mut state, items[index]) {
                                Ok(result) => {
                                    results.lock().expect("results lock")[index] = Some(result);
                                }
                                Err(e) => {
                                    failed.store(true, Ordering::Relaxed);
                                    return Err(e);
                                }
                            }
                        }
                    },
                )
            })
            .collect::<std::io::Result<Vec<_>>>()
            .context("failed to start a prefetch thread")?;
        for handle in handles {
            handle.join().unwrap_or_else(|_| Err(anyhow!("a prefetch thread panicked")))?;
        }
        Ok(())
    })?;
    results
        .into_inner()
        .expect("results lock")
        .into_iter()
        .map(|r| r.ok_or_else(|| anyhow!("a prefetch item was not processed")))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batch_alignments_in_a_crowded_block_are_not_far() {
        let mut ids: Vec<i64> = (2048..2048 + 100).collect();
        assert!(far_in_batch(&mut ids).is_empty());
    }

    #[test]
    fn batch_alignments_alone_in_their_block_are_far() {
        let mut ids: Vec<i64> = (2048..2048 + 100).collect();
        ids.extend([5_000_000, 9_000_123]);
        assert_eq!(far_in_batch(&mut ids), vec![5_000_000, 9_000_123]);
    }

    #[test]
    fn a_few_alignments_sharing_a_block_are_far_together() {
        let mut ids: Vec<i64> = (2048..2048 + 100).collect();
        ids.extend([5_000_000, 5_000_001, 5_000_010]);
        assert_eq!(far_in_batch(&mut ids), vec![5_000_000, 5_000_001, 5_000_010]);
    }

    #[test]
    fn unaligned_reads_are_never_far() {
        let mut ids = vec![0, 0, 0];
        assert!(far_in_batch(&mut ids).is_empty());
    }

    fn chunk(reads: &[(i64, &[u8])]) -> Chunk {
        let mut chunk = Chunk { ids: Vec::new(), starts: Vec::new(), bases: Vec::new() };
        for &(id, bases) in reads {
            chunk.ids.push(id);
            chunk.starts.push(chunk.bases.len());
            chunk.bases.extend_from_slice(bases);
        }
        chunk.starts.push(chunk.bases.len());
        chunk
    }

    #[test]
    fn prefetched_reads_are_found_by_id_across_chunks() {
        let reads = PrefetchedReads::from_chunks(vec![
            chunk(&[(3, b"ACGT"), (70, b"GG")]),
            chunk(&[(200, b"TTTAA")]),
        ]);
        assert_eq!(reads.get(3), Some(&b"ACGT"[..]));
        assert_eq!(reads.get(70), Some(&b"GG"[..]));
        assert_eq!(reads.get(200), Some(&b"TTTAA"[..]));
        assert_eq!(reads.len(), 3);
        assert_eq!(reads.total_bases(), 11);
    }

    #[test]
    fn ids_not_prefetched_are_absent() {
        let reads = PrefetchedReads::from_chunks(vec![chunk(&[(3, b"ACGT")])]);
        assert_eq!(reads.get(4), None);
        assert_eq!(reads.get(0), None);
        assert_eq!(reads.get(-1), None);
        assert_eq!(reads.get(1_000_000), None);
    }

    #[test]
    fn nothing_prefetched_finds_nothing() {
        let reads = PrefetchedReads::from_chunks(Vec::new());
        assert_eq!(reads.get(1), None);
        assert_eq!(reads.len(), 0);
    }

    #[test]
    fn parallel_chunks_keeps_item_order() {
        let items: Vec<usize> = (0..100).collect();
        let doubled = parallel_chunks(&items, vec![(); 4], |(), i| Ok(i * 2)).unwrap();
        assert_eq!(doubled, (0..100).map(|i| i * 2).collect::<Vec<_>>());
    }

    // A large stack frame is the point: 4 MiB is twice Rust's default for a thread, and well
    // under the 16 MiB given even if a debug build copies the array.
    #[allow(clippy::large_stack_arrays)]
    #[test]
    fn parallel_chunks_threads_have_room_for_deep_library_stacks() {
        let doubled = parallel_chunks(&[1, 2, 3], vec![(), ()], |(), item: i32| {
            let stack = std::hint::black_box([1u8; 4 << 20]);
            Ok(item * 2 * i32::from(stack[stack.len() - 1]))
        });
        assert_eq!(doubled.unwrap(), [2, 4, 6]);
    }

    #[test]
    fn parallel_chunks_reports_an_error() {
        let items: Vec<usize> = (0..100).collect();
        let result = parallel_chunks(&items, vec![(); 4], |(), i| {
            if i == 42 { Err(anyhow!("item 42 failed")) } else { Ok(i) }
        });
        assert!(result.unwrap_err().to_string().contains("item 42"));
    }
}
