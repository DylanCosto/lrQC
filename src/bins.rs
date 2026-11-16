//! Binning of read positions and lengths.

use std::sync::LazyLock;

/// Position bin widths: one base up to 100, then wider bins each decade, as FastQC
/// groups positions, extended to the lengths of long reads.
const POSITION_STEPS: [(u64, u64); 6] = [
    (100, 1),
    (1_000, 50),
    (10_000, 500),
    (100_000, 5_000),
    (1_000_000, 50_000),
    // Ultra-long reads: the longest reads ever sequenced are about 4 Mb.
    (10_000_000, 500_000),
];

/// Bins of read positions, counted from either end of a read.
#[derive(Debug, Clone)]
pub struct PositionBins {
    /// Start of each bin; the last bin is open-ended.
    starts: Vec<u64>,
}

/// The shared position bins.
pub static POSITIONS: LazyLock<PositionBins> = LazyLock::new(PositionBins::new);

impl PositionBins {
    fn new() -> Self {
        let mut starts = Vec::new();
        let mut at = 0;
        for (until, width) in POSITION_STEPS {
            while at < until {
                starts.push(at);
                at += width;
            }
        }
        starts.push(at);
        Self { starts }
    }

    /// Number of bins.
    pub fn len(&self) -> usize {
        self.starts.len()
    }

    /// The bin's first position (0-based) and, unless it is the open last bin, its end.
    pub fn range(&self, bin: usize) -> (u64, Option<u64>) {
        (self.starts[bin], self.starts.get(bin + 1).copied())
    }

    /// Calls `visit(bin, start, end)` for each bin overlapping positions `0..len`.
    pub fn for_each_segment(&self, len: usize, mut visit: impl FnMut(usize, usize, usize)) {
        for (bin, &start) in self.starts.iter().enumerate() {
            let start = usize::try_from(start).unwrap_or(usize::MAX);
            if start >= len {
                break;
            }
            let end = self.starts.get(bin + 1).map_or(len, |&next| {
                usize::try_from(next).unwrap_or(usize::MAX).min(len)
            });
            visit(bin, start, end);
        }
    }

    /// The bin holding `position`.
    pub fn bin_of(&self, position: u64) -> usize {
        self.starts.partition_point(|&s| s <= position) - 1
    }

    /// Bin width, taking the open last bin as one wide as the one before it.
    pub fn width(&self, bin: usize) -> u64 {
        match self.range(bin) {
            (start, Some(end)) => end - start,
            (start, None) => start - self.starts[bin - 1],
        }
    }
}

/// Log-scaled read length bins: this many per factor of ten.
pub const LENGTH_BINS_PER_DECADE: f64 = 20.0;

/// Number of log length bins, covering lengths up to 10 Mb.
pub const LENGTH_BINS: usize = 141;

/// The log bin of a read length.
pub fn length_bin(len: u64) -> usize {
    if len == 0 {
        return 0;
    }
    // Lengths are far below 2^52, so the conversion is exact.
    let bin = ((len as f64).log10() * LENGTH_BINS_PER_DECADE).floor();
    (bin.max(0.0) as usize).min(LENGTH_BINS - 1)
}

/// Lower edge of a log length bin.
pub fn length_bin_start(bin: usize) -> f64 {
    10f64.powf(bin as f64 / LENGTH_BINS_PER_DECADE)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn segments_cover_the_read_exactly_once() {
        for len in [0, 1, 99, 100, 101, 1234, 56_789, 2_000_000, 12_345_678] {
            let mut covered = 0;
            let mut last_end = 0;
            POSITIONS.for_each_segment(len, |bin, start, end| {
                assert_eq!(start, last_end, "gap before bin {bin} for length {len}");
                assert!(end > start);
                covered += end - start;
                last_end = end;
            });
            assert_eq!(covered, len);
        }
    }

    #[test]
    fn finds_bins() {
        assert_eq!(POSITIONS.bin_of(0), 0);
        assert_eq!(POSITIONS.bin_of(99), 99);
        assert_eq!(POSITIONS.bin_of(100), 100);
        assert_eq!(POSITIONS.bin_of(149), 100);
        assert_eq!(POSITIONS.bin_of(150), 101);
        assert_eq!(POSITIONS.bin_of(u64::MAX), POSITIONS.len() - 1);
        assert_eq!(POSITIONS.width(100), 50);
    }

    #[test]
    fn length_bins_are_monotonic() {
        assert_eq!(length_bin(1), 0);
        assert_eq!(length_bin(10), 20);
        assert_eq!(length_bin(10_000), 80);
        assert_eq!(length_bin(u64::MAX), LENGTH_BINS - 1);
        assert!(length_bin(1_000) < length_bin(1_200));
    }
}
