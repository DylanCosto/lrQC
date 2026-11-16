//! A reproducible sample of read ends, for duplication and overrepresented sequences.
//!
//! Each group keeps the reads with the smallest name hashes (a bottom-k sample). The sample
//! depends only on the read names, not on thread scheduling, and two samples merge exactly.

use std::collections::BinaryHeap;

use crate::read::BASE_N;

/// Reads sampled per group; FastQC samples 100,000 short reads, but long-read groups are
/// many and one percent of 10,000 is still 100 reads.
pub const SAMPLE_SIZE: usize = 10_000;

/// Bases taken from each end of a read, as FastQC truncates reads longer than 75 bases to 50.
pub const END_LEN: usize = 50;

/// A read's two ends, 2-bit packed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Ends {
    pub start: u128,
    pub end: u128,
}

/// Bottom-k sample of read ends.
#[derive(Debug, Clone, Default)]
pub struct EndSample {
    heap: BinaryHeap<(u64, Ends)>,
}

impl EndSample {
    /// Offers a read; it is kept if its name hash is among the smallest seen.
    pub fn offer(&mut self, name_hash: u64, bases: &[u8]) {
        if bases.len() < END_LEN {
            return;
        }
        if self.heap.len() >= SAMPLE_SIZE
            && self.heap.peek().is_some_and(|&(top, _)| name_hash >= top)
        {
            return;
        }
        let (Some(start), Some(end)) = (
            pack(&bases[..END_LEN]),
            pack(&bases[bases.len() - END_LEN..]),
        ) else {
            return;
        };
        self.heap.push((name_hash, Ends { start, end }));
        if self.heap.len() > SAMPLE_SIZE {
            self.heap.pop();
        }
    }

    /// Combines two samples into the sample of their union.
    pub fn merge(&mut self, other: EndSample) {
        for item in other.heap {
            self.heap.push(item);
            if self.heap.len() > SAMPLE_SIZE {
                self.heap.pop();
            }
        }
    }

    /// The sampled reads' ends.
    pub fn ends(&self) -> impl Iterator<Item = Ends> + '_ {
        self.heap.iter().map(|&(_, ends)| ends)
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }
}

/// 2-bit packs up to 64 bases, or `None` if any is N.
fn pack(bases: &[u8]) -> Option<u128> {
    bases.iter().try_fold(0u128, |acc, &c| {
        (c != BASE_N).then_some((acc << 2) | u128::from(c))
    })
}

/// Unpacks [`END_LEN`] bases packed by [`EndSample::offer`].
pub fn unpack(packed: u128) -> String {
    (0..END_LEN)
        .rev()
        .map(|i| char::from(b"ACGT"[((packed >> (2 * i)) & 3) as usize]))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_smallest_hashes_regardless_of_order() {
        let bases = vec![1u8; 60];
        let mut forward = EndSample::default();
        let mut split_a = EndSample::default();
        let mut split_b = EndSample::default();
        for hash in 0..(3 * SAMPLE_SIZE as u64) {
            let scrambled = hash.wrapping_mul(0x9e37_79b9_7f4a_7c15);
            forward.offer(scrambled, &bases);
            if hash % 2 == 0 {
                split_a.offer(scrambled, &bases);
            } else {
                split_b.offer(scrambled, &bases);
            }
        }
        split_a.merge(split_b);
        let mut a: Vec<u64> = forward.heap.iter().map(|x| x.0).collect();
        let mut b: Vec<u64> = split_a.heap.iter().map(|x| x.0).collect();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
        assert_eq!(a.len(), SAMPLE_SIZE);
    }

    #[test]
    fn round_trips_ends() {
        let seq = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAAC";
        let codes: Vec<u8> = seq
            .bytes()
            .map(|b| crate::read::ASCII_TO_CODE[usize::from(b)])
            .collect();
        assert_eq!(unpack(pack(&codes).expect("no N")), seq);
    }
}
