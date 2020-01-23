use crate::ivstore::IntervalLike;
use std::cmp::Ordering::{self};

/// Represent a range from [start, stop)
/// Inclusive start, exclusive of stop
#[derive(Debug)]
pub struct Interval<T> {
    pub start: u32,
    pub stop: u32,
    pub val: T,
}

impl<T> IntervalLike<T> for Interval<T> {
    /// Compute the intsect between two intervals
    #[inline]
    fn intersect(&self, other: &Interval<T>) -> u32 {
        std::cmp::min(self.stop, other.stop)
            .checked_sub(std::cmp::max(self.start, other.start))
            .unwrap_or(0)
    }

    /// Check if two intervals overlap
    #[inline]
    fn overlap(&self, start: u32, stop: u32) -> bool {
        self.start < stop && self.stop > start
    }

    #[inline]
    fn start(&self) -> u32 {
        self.start
    }

    #[inline]
    fn stop(&self) -> u32 {
        self.stop
    }

    #[inline]
    fn val(&self) -> &T {
        &self.val
    }

    #[inline]
    fn set_start(&mut self, new: u32) {
        self.start = new;
    }

    #[inline]
    fn set_stop(&mut self, new: u32) {
        self.stop = new;
    }

    #[inline]
    fn set_val(&mut self, new: T) {
        self.val = new;
    }
}

impl<T> Ord for Interval<T> {
    #[inline]
    fn cmp(&self, other: &Interval<T>) -> Ordering {
        if self.start < other.start {
            Ordering::Less
        } else if other.start < self.start {
            Ordering::Greater
        } else {
            self.stop.cmp(&other.stop)
        }
    }
}
impl<T> Eq for Interval<T> {}

impl<T> PartialOrd for Interval<T> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl<T> PartialEq for Interval<T> {
    #[inline]
    fn eq(&self, other: &Interval<T>) -> bool {
        self.start == other.start && self.stop == other.stop
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod test {
    use super::*;
    type Iv = Interval<u32>;
    #[test]
    fn test_interval_intersects() {
        let i1 = Iv{start: 70, stop: 120, val: 0}; // max_len = 50
        let i2 = Iv{start: 10, stop: 15, val: 0};
        let i3 = Iv{start: 10, stop: 15, val: 0}; // exact overlap
        let i4 = Iv{start: 12, stop: 15, val: 0}; // inner overlap
        let i5 = Iv{start: 14, stop: 16, val: 0}; // overlap end
        let i6 = Iv{start: 40, stop: 50, val: 0};
        let i7 = Iv{start: 50, stop: 55, val: 0};
        let i_8 = Iv{start: 60, stop: 65, val: 0};
        let i9 = Iv{start: 68, stop: 71, val: 0}; // overlap start
        let i10 = Iv{start: 70, stop: 75, val: 0};

        assert_eq!(i2.intersect(&i3), 5); // exact match
        assert_eq!(i2.intersect(&i4), 3); // inner intersect
        assert_eq!(i2.intersect(&i5), 1); // end intersect
        assert_eq!(i9.intersect(&i10), 1); // start intersect
        assert_eq!(i7.intersect(&i_8), 0); // no intersect
        assert_eq!(i6.intersect(&i7), 0); // no intersect stop = start
        assert_eq!(i1.intersect(&i10), 5); // inner intersect at start
    }
}
