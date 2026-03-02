//! FxHashMap-based mate-pair cache.
//!
//! Stores mate information keyed by alignment ID for resolving mate fields
//! (RNEXT, PNEXT, TLEN) in paired-end reads. Cleared between references
//! in default (non-region) mode.

use rustc_hash::FxHashMap;

/// Information stored about an alignment for its mate to look up.
///
/// The mate cache is cleared between references, so all cached mates are
/// on the same reference as the current read — `ref_name` is not stored.
#[derive(Debug, Clone)]
pub struct MateInfo {
    /// 0-based reference position of this alignment.
    pub ref_pos: i32,
    /// Template length.
    pub tlen: i32,
}

/// Cache for resolving mate-pair information in paired-end reads.
///
/// When we encounter an alignment, we store its info keyed by its own
/// alignment row ID. When we later encounter its mate (which has
/// `MATE_ALIGN_ID` pointing back to this row), we can look up the
/// mate's reference name, position, etc.
pub struct MateCache {
    map: FxHashMap<i64, MateInfo>,
}

impl MateCache {
    /// Create a new empty mate cache.
    pub fn new() -> Self {
        Self { map: FxHashMap::default() }
    }

    /// Insert mate info for a given alignment row ID.
    pub fn insert(&mut self, align_id: i64, info: MateInfo) {
        self.map.insert(align_id, info);
    }

    /// Look up and remove mate info for a given alignment row ID.
    ///
    /// Uses remove (not get) because each mate pair is only looked up once.
    pub fn take(&mut self, mate_align_id: i64) -> Option<MateInfo> {
        self.map.remove(&mate_align_id)
    }

    /// Clear the cache (typically between references).
    pub fn clear(&mut self) {
        self.map.clear();
    }

    /// Number of entries in the cache.
    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.map.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insert_and_take() {
        let mut cache = MateCache::new();
        let info = MateInfo { ref_pos: 1000, tlen: 300 };
        cache.insert(42, info);

        let mate = cache.take(42).expect("should find mate");
        assert_eq!(mate.ref_pos, 1000);
        assert_eq!(mate.tlen, 300);

        // Second take should return None (removed).
        assert!(cache.take(42).is_none());
    }

    #[test]
    fn test_take_missing() {
        let mut cache = MateCache::new();
        assert!(cache.take(999).is_none());
    }

    #[test]
    fn test_clear() {
        let mut cache = MateCache::new();
        cache.insert(1, MateInfo { ref_pos: 0, tlen: 0 });
        cache.insert(2, MateInfo { ref_pos: 0, tlen: 0 });
        assert_eq!(cache.len(), 2);

        cache.clear();
        assert_eq!(cache.len(), 0);
        assert!(cache.take(1).is_none());
    }
}
