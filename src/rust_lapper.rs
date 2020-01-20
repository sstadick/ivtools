//! This module provides a simple data-structure for fast interval searches.
//! ## Features
//! - Extremely fast on most genomic datasets. (3-4x faster than other methods)
//! - Extremely fast on in order queries. (10x faster than other methods)
//! - Extremely fast intersections count method based on the
//! [BITS](https://arxiv.org/pdf/1208.3407.pdf) algorithm
//! - Parallel friendly. Queries are on an immutable structure, even for seek
//! - Consumer / Adapter paradigm, Iterators are returned and serve as the main API for interacting
//! with the lapper
//!
//! ## Details:
//!
//! ```text
//!       0  1  2  3  4  5  6  7  8  9  10 11
//! (0,10]X  X  X  X  X  X  X  X  X  X
//! (2,5]       X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (5,9]                X  X  X  X
//! (8,11]                        X  X  X
//!
//! Query: (8, 11]
//! Answer: ((0,10], (5,9], (8,11])
//! ```
//!
//! Most interaction with this crate will be through the [`Lapper``](struct.Lapper.html)` struct
//! The main methods are [`find`](struct.Lapper.html#method.find),
//! [`seek`](struct.Lapper.html#method.seek), and [`count`](struct.Lapper.html#method.count)
//! where both `seek` and `count` are special cases allowing for very fast queries in certain scenarios.
//!
//! The overlap function for this assumes a zero based genomic coordinate system. So [start, stop)
//! is not inclusive of the stop position for neither the queries, nor the Intervals.
//!
//! Lapper does not use an interval tree, instead, it operates on the assumtion that most intervals are
//! of similar length; or, more exactly, that the longest interval in the set is not long compred to
//! the average distance between intervals.
//!
//! For cases where this holds true (as it often does with genomic data), we can sort by start and
//! use binary search on the starts, accounting for the length of the longest interval. The advantage
//! of this approach is simplicity of implementation and speed. In realistic tests queries returning
//! the overlapping intervals are 1000 times faster than brute force and queries that merely check
//! for the overlaps are > 5000 times faster.
//!
//! When this is not the case, if possible in your scenario, use merge_overlaps first, and then use
//! `find` or `seek`. The `count` method will be fast in all scenarios.
//!
//! # Examples
//!
//! ```rust
//!    use rust_lapper::{Interval, Lapper};
//!    use std::cmp;
//!    type Iv = Interval<u32>;
//!
//!    // create some fake data
//!    let data: Vec<Iv> = (0..20).step_by(5).map(|x| Iv{start: x, stop: x + 2, val: 0}).collect();
//!    println!("{:#?}", data);
//!
//!    // make lapper structure
//!    let laps = Lapper::new(data);
//!
//!    assert_eq!(laps.find(6, 11).next(), Some(&Iv{start: 5, stop: 7, val: 0}));
//!    
//!    // Demonstration of seek function. By passing in the &mut cursor, seek can have thread local
//!    // cursors going
//!    let mut sim: u32= 0;
//!    let mut cursor = 0;
//!    // Calculate the overlap between the query and the found intervals, sum total overlap
//!    for i in (0..10).step_by(3) {
//!        sim += laps
//!            .seek(i, i + 2, &mut cursor)
//!            .map(|iv| cmp::min(i + 2, iv.stop) - cmp::max(i, iv.start))
//!            .sum::<u32>();
//!    }
//!    assert_eq!(sim, 4);
//! ```
use crate::interval::Interval;
use crate::ivstore::IvStore;
use std::collections::VecDeque;

/// Primary object of the library. The public intervals holds all the intervals and can be used for
/// iterating / pulling values out of the tree.
#[derive(Debug)]
pub struct Lapper<T> {
    /// List of intervasl
    pub intervals: Vec<Interval<T>>,
    /// The length of the longest interval
    max_len: u32,
    /// A cursor to hold the position in the list in between searches with `seek` method
    cursor: usize,
}

impl<T> Lapper<T> {
    /// Determine the first index that we should start checking for overlaps for via a binary
    /// search.
    #[inline]
    pub fn lower_bound(start: u32, intervals: &[Interval<T>]) -> usize {
        let mut size = intervals.len();
        let mut low = 0;

        while size > 0 {
            let half = size / 2;
            let other_half = size - half;
            let probe = low + half;
            let other_low = low + other_half;
            let v = &intervals[probe];
            size = half;
            low = if v.start < start { other_low } else { low }
        }
        low
    }

    /// Find all intevals that overlap start .. stop. This method will work when queries
    /// to this lapper are in sorted (start) order. It uses a linear search from the last query
    /// instead of a binary search. A reference to a cursor must be passed in. This reference will
    /// be modified and should be reused in the next query. This allows seek to not need to make
    /// the lapper object mutable, and thus use the same lapper accross threads.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<bool>>>());
    /// let mut cursor = 0;
    /// for i in lapper.iter() {
    ///    assert_eq!(lapper.seek(i.start, i.stop, &mut cursor).count(), 1);
    /// }
    /// ```
    #[inline]
    pub fn seek<'a>(&'a self, start: u32, stop: u32, cursor: &mut usize) -> IterFind<'a, T> {
        if *cursor == 0 || (*cursor < self.intervals.len() && self.intervals[*cursor].start > start)
        {
            *cursor = Self::lower_bound(
                start.checked_sub(self.max_len).unwrap_or(0),
                &self.intervals,
            );
        }

        while *cursor + 1 < self.intervals.len()
            && self.intervals[*cursor + 1].start < start.checked_sub(self.max_len).unwrap_or(0)
        {
            *cursor += 1;
        }

        IterFind {
            inner: self,
            off: *cursor,
            end: self.intervals.len(),
            start,
            stop,
        }
    }
}

impl<'a, T> IvStore<'a, T> for Lapper<T>
where
    T: 'a,
{
    type IvSearchIterator = IterFind<'a, T>;
    type IvIterator = IterLapper<'a, T>;
    type SelfU32 = Lapper<u32>;

    /// Create a new instance of Lapper by passing in a vector of Intervals. This vector will
    /// immediately be sorted by start order.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, stop: x + 10, val: true})
    ///                   .collect::<Vec<Interval<bool>>>();
    /// let lapper = Lapper::new(data);
    /// ```
    fn new(mut intervals: Vec<Interval<T>>) -> Self {
        intervals.sort();
        let mut max_len = 0;
        for interval in intervals.iter() {
            let i_len = interval.stop.checked_sub(interval.start).unwrap_or(0);
            if i_len > max_len {
                max_len = i_len;
            }
        }
        Lapper {
            intervals,
            max_len,
            cursor: 0,
        }
    }

    /// Get the number over intervals in Lapper
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, stop: x + 10, val: true})
    ///                   .collect::<Vec<Interval<bool>>>();
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.len(), 4);
    /// ```
    #[inline]
    fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Check if lapper is empty
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data: Vec<Interval<bool>> = vec![];
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.is_empty(), true);
    /// ```
    #[inline]
    fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Return an iterator over the intervals in Lapper
    #[inline]
    fn iter(&self) -> IterLapper<T> {
        IterLapper {
            inner: self,
            pos: 0,
        }
    }

    /// Merge any intervals that overlap with eachother within the Lapper. This is an easy way to
    /// speed up queries. It returns a new lapper
    fn merge_overlaps(mut self) -> Lapper<u32> {
        let mut stack: VecDeque<&mut Interval<T>> = VecDeque::new();
        let mut ivs = self.intervals.iter_mut();
        let intervals: Vec<Interval<u32>> = if let Some(first) = ivs.next() {
            stack.push_back(first);
            for interval in ivs {
                let mut top = stack.pop_back().unwrap();
                if top.stop < interval.start {
                    stack.push_back(top);
                    stack.push_back(interval);
                } else if top.stop < interval.stop {
                    top.stop = interval.stop;
                    //stack.pop_back();
                    stack.push_back(top);
                } else {
                    // they were equal
                    stack.push_back(top);
                }
            }
            stack
                .into_iter()
                .map(|x| Interval {
                    start: x.start,
                    stop: x.stop,
                    // val: x.val.clone(),
                    val: 0,
                })
                .collect()
        } else {
            vec![]
        };
        Lapper::new(intervals)
        // Fix the starts and stops used by counts
        // let (mut starts, mut stops): (Vec<_>, Vec<_>) =
        //     intervals.iter().map(|x| (x.start, x.stop)).unzip();
        // starts.sort();
        // stops.sort();
        // self.starts = starts;
        // self.stops = stops;
    }

    /// Find all intervals that overlap start .. stop
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<bool>>>());
    /// assert_eq!(lapper.find(5, 11).count(), 2);
    /// ```
    #[inline]
    fn find(&self, start: u32, stop: u32) -> IterFind<T> {
        IterFind {
            inner: self,
            off: Self::lower_bound(
                start.checked_sub(self.max_len).unwrap_or(0),
                &self.intervals,
            ),
            end: self.intervals.len(),
            start,
            stop,
        }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, T>
where
    T: 'a,
{
    inner: &'a Lapper<T>,
    off: usize,
    end: usize,
    start: u32,
    stop: u32,
}

impl<'a, T> Iterator for IterFind<'a, T> {
    type Item = &'a Interval<T>;

    #[inline]
    // interval.start < stop && interval.stop > start
    fn next(&mut self) -> Option<Self::Item> {
        while self.off < self.inner.intervals.len() {
            //let mut generator = self.inner.intervals[self.off..].iter();
            //while let Some(interval) = generator.next() {
            let interval = &self.inner.intervals[self.off];
            self.off += 1;
            if interval.overlap(self.start, self.stop) {
                return Some(interval);
            } else if interval.start >= self.stop {
                break;
            }
        }
        None
    }
}

/// Lapper Iterator
pub struct IterLapper<'a, T>
where
    T: 'a,
{
    inner: &'a Lapper<T>,
    pos: usize,
}

impl<'a, T> Iterator for IterLapper<'a, T> {
    type Item = &'a Interval<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.intervals.len() {
            None
        } else {
            self.pos += 1;
            self.inner.intervals.get(self.pos - 1)
        }
    }
}

impl<T> IntoIterator for Lapper<T> {
    type Item = Interval<T>;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intervals.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Lapper<T> {
    type Item = &'a Interval<T>;
    type IntoIter = std::slice::Iter<'a, Interval<T>>;

    fn into_iter(self) -> std::slice::Iter<'a, Interval<T>> {
        self.intervals.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut Lapper<T> {
    type Item = &'a mut Interval<T>;
    type IntoIter = std::slice::IterMut<'a, Interval<T>>;

    fn into_iter(self) -> std::slice::IterMut<'a, Interval<T>> {
        self.intervals.iter_mut()
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    type Iv = Interval<u32>;
    fn setup_nonoverlapping() -> Lapper<u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(20)
            .map(|x| Iv {
                start: x,
                stop: x + 10,
                val: 0,
            })
            .collect();
        let lapper = Lapper::new(data);
        lapper
    }
    fn setup_overlapping() -> Lapper<u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(10)
            .map(|x| Iv {
                start: x,
                stop: x + 15,
                val: 0,
            })
            .collect();
        let lapper = Lapper::new(data);
        lapper
    }
    fn setup_badlapper() -> Lapper<u32> {
        let data: Vec<Iv> = vec![
            Iv{start: 70, stop: 120, val: 0}, // max_len = 50
            Iv{start: 10, stop: 15, val: 0},
            Iv{start: 10, stop: 15, val: 0}, // exact overlap
            Iv{start: 12, stop: 15, val: 0}, // inner overlap
            Iv{start: 14, stop: 16, val: 0}, // overlap end
            Iv{start: 40, stop: 45, val: 0},
            Iv{start: 50, stop: 55, val: 0},
            Iv{start: 60, stop: 65, val: 0},
            Iv{start: 68, stop: 71, val: 0}, // overlap start
            Iv{start: 70, stop: 75, val: 0},
        ];
        let lapper = Lapper::new(data);
        lapper
    }
    fn setup_single() -> Lapper<u32> {
        let data: Vec<Iv> = vec![Iv {
            start: 10,
            stop: 35,
            val: 0,
        }];
        let lapper = Lapper::new(data);
        lapper
    }

    // Test that a query stop that hits an interval start returns no interval
    #[test]
    fn test_query_stop_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(15, 20).next());
        assert_eq!(None, lapper.seek(15, 20, &mut cursor).next());
    }

    // Test that a query start that hits an interval end returns no interval
    #[test]
    fn test_query_start_interval_stop() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(30, 35).next());
        assert_eq!(None, lapper.seek(30, 35, &mut cursor).next());
    }

    // Test that a query that overlaps the start of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 25).next());
        assert_eq!(Some(&expected), lapper.seek(15, 25, &mut cursor).next());
    }

    // Test that a query that overlaps the stop of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_stop() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(25, 35).next());
        assert_eq!(Some(&expected), lapper.seek(25, 35, &mut cursor).next());
    }

    // Test that a query that is enveloped by interval returns interval
    #[test]
    fn test_interval_envelops_query() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(22, 27).next());
        assert_eq!(Some(&expected), lapper.seek(22, 27, &mut cursor).next());
    }

    // Test that a query that envolops an interval returns that interval
    #[test]
    fn test_query_envolops_interval() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            stop: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 35).next());
        assert_eq!(Some(&expected), lapper.seek(15, 35, &mut cursor).next());
    }

    #[test]
    fn test_overlapping_intervals() {
        let lapper = setup_overlapping();
        let mut cursor = 0;
        let e1 = Iv {
            start: 0,
            stop: 15,
            val: 0,
        };
        let e2 = Iv {
            start: 10,
            stop: 25,
            val: 0,
        };
        assert_eq!(vec![&e1, &e2], lapper.find(8, 20).collect::<Vec<&Iv>>());
        assert_eq!(
            vec![&e1, &e2],
            lapper.seek(8, 20, &mut cursor).collect::<Vec<&Iv>>()
        );
    }

    #[test]
    fn test_merge_overlaps() {
        let lapper = setup_badlapper();
        let expected: Vec<&Iv> = vec![
            &Iv{start: 10, stop: 16, val: 0},
            &Iv{start: 40, stop: 45, val: 0},
            &Iv{start: 50, stop: 55, val: 0},
            &Iv{start: 60, stop: 65, val: 0},
            &Iv{start: 68, stop: 120, val: 0}, // max_len = 50
        ];
        let new_lapper = lapper.merge_overlaps();
        assert_eq!(expected, new_lapper.iter().collect::<Vec<&Iv>>());
        
    }

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


    #[test]
    fn test_find_overlaps_in_large_intervals() {
        let data1: Vec<Iv> = vec![
            Iv{start: 0, stop: 8, val: 0},
            Iv{start: 1, stop: 10, val: 0}, 
            Iv{start: 2, stop: 5, val: 0}, 
            Iv{start: 3, stop: 8, val: 0},
            Iv{start: 4, stop: 7, val: 0},
            Iv{start: 5, stop: 8, val: 0},
            Iv{start: 8, stop: 8, val: 0},
            Iv{start: 9, stop: 11, val: 0},
            Iv{start: 10, stop: 13, val: 0},
            Iv{start: 100, stop: 200, val: 0},
            Iv{start: 110, stop: 120, val: 0},
            Iv{start: 110, stop: 124, val: 0},
            Iv{start: 111, stop: 160, val: 0},
            Iv{start: 150, stop: 200, val: 0},
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.find(8, 11).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 1, stop: 10, val: 0}, 
            &Iv{start: 9, stop: 11, val: 0},
            &Iv{start: 10, stop: 13, val: 0},
        ]);
        let found = lapper.find(145, 151).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 100, stop: 200, val: 0},
            &Iv{start: 111, stop: 160, val: 0},
            &Iv{start: 150, stop: 200, val: 0},
        ]);

    }

  

    
    // BUG TESTS - these are tests that came from real life

    // Test that it's not possible to induce index out of bounds by pushing the cursor past the end
    // of the lapper.
    #[test]
    fn test_seek_over_len() {
        let lapper = setup_nonoverlapping();
        let single = setup_single();
        let mut cursor: usize = 0;

        for interval in lapper.iter() {
            for o_interval in single.seek(interval.start, interval.stop, &mut cursor) {
                println!("{:#?}", o_interval);
            }
        }
    }

    // Test that if lower_bound puts us before the first match, we still return a match
    #[test]
    fn test_find_over_behind_first_match() {
        let lapper = setup_badlapper();
        let e1 = Iv {start: 50, stop: 55, val: 0};
        let found = lapper.find(50, 55).next();
        assert_eq!(found, Some(&e1));
    }

    // Test that seeking for all intervals over self == len(self)
    #[test]
    fn test_seek_over_nonoverlapping() {
        let lapper = setup_nonoverlapping();
        let mut total = 0;
        let mut cursor: usize = 0;
        for iv in lapper.iter() {
            for _ in lapper.seek(iv.start, iv.stop, &mut cursor) {
                total += 1;
            }
        }
        assert_eq!(lapper.len(), total);
    }

    // Test that finding for all intervals over self == len(self)
    #[test]
    fn test_find_over_nonoverlapping() {
        let lapper = setup_nonoverlapping();
        let mut total = 0;
        for iv in lapper.iter() {
            for _ in lapper.find(iv.start, iv.stop) {
                total += 1;
            }
        }
        assert_eq!(lapper.len(), total);
    }

    // When there is a very long interval that spans many little intervals, test that the little
    // intevals still get returne properly
    #[test]
    fn test_bad_skips() {
        let data = vec![
            Iv{start:25264912, stop: 25264986, val: 0},	
            Iv{start:27273024, stop: 27273065	, val: 0},
            Iv{start:27440273, stop: 27440318	, val: 0},
            Iv{start:27488033, stop: 27488125	, val: 0},
            Iv{start:27938410, stop: 27938470	, val: 0},
            Iv{start:27959118, stop: 27959171	, val: 0},
            Iv{start:28866309, stop: 33141404	, val: 0},
        ];
        let lapper = Lapper::new(data);

        let found = lapper.find(28974798, 33141355).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start:28866309, stop: 33141404	, val: 0},
        ]);
    }
}
