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
//!    use ivtools::{Interval, Lapper, IvStore};
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
use crate::ivstore::{IntervalLike, IvStore};
use std::collections::VecDeque;
use std::marker::PhantomData;

/// Primary object of the library. The public intervals holds all the intervals and can be used for
/// iterating / pulling values out of the tree.
#[derive(Debug)]
pub struct Lapper<T, I>
where
    T: Default,
    I: IntervalLike<T>,
{
    /// List of intervasl
    pub intervals: Vec<I>,
    /// The length of the longest interval
    max_len: u32,
    /// A cursor to hold the position in the list in between searches with `seek` method
    /// The div number created by the universe that is the length of my intervals (see van emde
    /// boas)
    universe_div: usize,
    /// An index of clusters of size sqrt(intervals.len()) as usize
    index: Vec<Cluster>,
    cursor: usize,
    phantom_type: PhantomData<T>,
    // total_vecs: usize,
    // starts: Vec<u32x8>,
    // stops: Vec<u32x8>,
    // max_lens: Vec<u32>,
    // max_ends: Vec<u32>,
}

/// Index of my intervals, kind of like a van emde boas concept of a cluster
#[derive(Debug)]
struct Cluster {
    // The smallest start value of any interval in the cluster
    min_start: u32,
    // The largest stop value of any interval in the cluster
    max_stop: u32,
    // The index in intervals that the cluster starts at
    start_index: usize,
    // The index in intervals that the cluster stop at
    stop_index: usize,
    // The max length of any interval in the cluster
    max_len: u32,
}

impl<T, I> Lapper<T, I>
where
    T: Default,
    I: IntervalLike<T>,
{
    /// Determine the first index that we should start checking for overlaps for via a binary
    /// search.
    #[inline]
    pub fn lower_bound(start: u32, intervals: &[I]) -> usize {
        let mut size = intervals.len();
        let mut low = 0;

        while size > 0 {
            let half = size / 2;
            let other_half = size - half;
            let probe = low + half;
            let other_low = low + other_half;
            let v = &intervals[probe];
            size = half;
            low = if v.start() < start { other_low } else { low }
        }
        low
    }

    pub fn super_find<'a>(&'a self, start: u32, stop: u32) -> Vec<&'a I> {
        let mut result = vec![];
        'index: for cluster in self.index.iter() {
            // let cluster = &self.index[cluster_offset];
            // cluster_offset += 1;
            if start < cluster.max_stop && stop > cluster.min_start {
                // there is at least one match in the cluster
                let mut offset = Self::lower_bound(
                    start.checked_sub(cluster.max_len).unwrap_or(0),
                    &self.intervals[cluster.start_index..cluster.stop_index],
                ) + cluster.start_index;
                while offset < cluster.stop_index && start < cluster.max_stop {
                    // missing an early exit criteria here
                    let interval = &self.intervals[offset];
                    offset += 1;
                    if interval.overlap(start, stop) {
                        result.push(interval);
                    } else if interval.start() >= stop {
                        break 'index;
                    }
                }
            }
        }
        result
    }

    /// Create the index based on a set of intervals
    fn build_index(intervals: &Vec<I>) -> (usize, Vec<Cluster>) {
        let interval_len = intervals.len();
        let clusters = (interval_len as f64).log2().floor() as usize + 1;
        let mut index: Vec<Cluster> = Vec::with_capacity(interval_len / clusters + 1);
        let mut counter = 0;
        for chunk in intervals.chunks(clusters) {
            let first = chunk.first().unwrap();
            let local_low_start = first.start();
            let local_high_stop = chunk
                .iter()
                .fold(first.stop(), |high, x| std::cmp::max(high, x.stop()));
            let local_max_len = chunk
                .iter()
                .fold(first.stop() - first.start(), |max_len, x| {
                    std::cmp::max(max_len, x.stop() - x.start())
                });
            index.push(Cluster {
                min_start: local_low_start,
                max_stop: local_high_stop,
                start_index: clusters * counter,
                stop_index: clusters * counter + chunk.len(),
                max_len: local_max_len,
            });
            counter += 1;
        }
        ((interval_len as f64).sqrt().ceil() as usize, index)
    }

    /// Find all intevals that overlap start .. stop. This method will work when queries
    /// to this lapper are in sorted (start) order. It uses a linear search from the last query
    /// instead of a binary search. A reference to a cursor must be passed in. This reference will
    /// be modified and should be reused in the next query. This allows seek to not need to make
    /// the lapper object mutable, and thus use the same lapper accross threads.
    /// ```
    /// use ivtools::{Lapper, Interval, IvStore};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<bool>>>());
    /// let mut cursor = 0;
    /// for i in lapper.iter() {
    ///    assert_eq!(lapper.seek(i.start, i.stop, &mut cursor).count(), 1);
    /// }
    /// ```
    #[inline]
    pub fn seek<'a>(&'a self, start: u32, stop: u32, cursor: &mut usize) -> IterFind<'a, T, I> {
        if *cursor == 0
            || (*cursor < self.intervals.len() && self.intervals[*cursor].start() > start)
        {
            *cursor = Self::lower_bound(
                start.checked_sub(self.max_len).unwrap_or(0),
                &self.intervals,
            );
        }

        while *cursor + 1 < self.intervals.len()
            && self.intervals[*cursor + 1].start() < start.checked_sub(self.max_len).unwrap_or(0)
        {
            *cursor += 1;
        }

        IterFind {
            inner: self,
            off: *cursor,
            end: self.intervals.len(),
            start,
            stop,
            phantom_type: PhantomData,
        }
    }
}

impl<'a, T, I> IvStore<'a, T, I> for Lapper<T, I>
where
    T: 'a + Default,
    I: IntervalLike<T> + 'a,
{
    type IvSearchIterator = IterFind<'a, T, I>;
    type IvIterator = IterLapper<'a, T, I>;

    /// Create a new instance of Lapper by passing in a vector of Intervals. This vector will
    /// immediately be sorted by start order.
    /// ```
    /// use ivtools::{Lapper, Interval, IvStore};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, stop: x + 10, val: true})
    ///                   .collect::<Vec<Interval<bool>>>();
    /// let lapper = Lapper::new(data);
    /// ```
    fn new(mut intervals: Vec<I>) -> Self {
        intervals.sort();
        let mut max_len = 0;
        for interval in intervals.iter() {
            let i_len = interval.stop().checked_sub(interval.start()).unwrap_or(0);
            if i_len > max_len {
                max_len = i_len;
            }
        }
        // create index
        let (universe_div, index) = Self::build_index(&intervals);

        Lapper {
            intervals,
            max_len,
            cursor: 0,
            phantom_type: PhantomData,
            universe_div: universe_div,
            index: index,
        }
    }

    /// Get the number over intervals in Lapper
    /// ```
    /// use ivtools::{Lapper, Interval, IvStore};
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
    /// use ivtools::{Lapper, Interval, IvStore};
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
    fn iter(&self) -> IterLapper<T, I> {
        IterLapper {
            inner: self,
            pos: 0,
            phantom_type: PhantomData,
        }
    }

    /// Merge any intervals that overlap with eachother within the Lapper. This is an easy way to
    /// speed up queries. It returns a new lapper
    fn merge_overlaps(mut self) -> Self {
        let mut stack: VecDeque<&mut I> = VecDeque::new();
        let mut ivs = self.intervals.iter_mut();
        let intervals = if let Some(first) = ivs.next() {
            stack.push_back(first);
            for interval in ivs {
                let top = stack.pop_back().unwrap();
                if top.stop() < interval.start() {
                    stack.push_back(top);
                    stack.push_back(interval);
                } else if top.stop() < interval.stop() {
                    top.set_stop(interval.stop());
                    //stack.pop_back();
                    stack.push_back(top);
                } else {
                    // they were equal
                    stack.push_back(top);
                }
            }
            stack
                .into_iter()
                .map(|x| I::new(x.start(), x.stop(), Default::default()))
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
    /// use ivtools::{Lapper, Interval, IvStore};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<bool>>>());
    /// assert_eq!(lapper.find(5, 11).count(), 2);
    /// ```
    #[inline]
    fn find(&self, start: u32, stop: u32) -> IterFind<T, I> {
        IterFind {
            inner: self,
            off: Self::lower_bound(
                start.checked_sub(self.max_len).unwrap_or(0),
                &self.intervals,
            ),
            end: self.intervals.len(),
            start,
            stop,
            phantom_type: PhantomData,
        }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, T, I>
where
    T: 'a + Default,
    I: IntervalLike<T>,
{
    inner: &'a Lapper<T, I>,
    off: usize,
    end: usize,
    start: u32,
    stop: u32,
    phantom_type: PhantomData<T>,
}

impl<'a, T, I> Iterator for IterFind<'a, T, I>
where
    T: Default,
    I: IntervalLike<T>,
{
    type Item = &'a I;

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
            } else if interval.start() >= self.stop {
                break;
            }
        }
        None
    }
}

/// Lapper Iterator
pub struct IterLapper<'a, T, I>
where
    T: 'a + Default,
    I: IntervalLike<T>,
{
    inner: &'a Lapper<T, I>,
    pos: usize,
    phantom_type: PhantomData<T>,
}

impl<'a, T: Default, I: IntervalLike<T>> Iterator for IterLapper<'a, T, I> {
    type Item = &'a I;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.intervals.len() {
            None
        } else {
            self.pos += 1;
            self.inner.intervals.get(self.pos - 1)
        }
    }
}

impl<T, I> IntoIterator for Lapper<T, I>
where
    T: Default,
    I: IntervalLike<T>,
{
    type Item = I;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intervals.into_iter()
    }
}

impl<'a, T, I> IntoIterator for &'a Lapper<T, I>
where
    T: Default,
    I: IntervalLike<T>,
{
    type Item = &'a I;
    type IntoIter = std::slice::Iter<'a, I>;

    fn into_iter(self) -> std::slice::Iter<'a, I> {
        self.intervals.iter()
    }
}

impl<'a, T, I> IntoIterator for &'a mut Lapper<T, I>
where
    T: Default,
    I: IntervalLike<T>,
{
    type Item = &'a mut I;
    type IntoIter = std::slice::IterMut<'a, I>;

    fn into_iter(self) -> std::slice::IterMut<'a, I> {
        self.intervals.iter_mut()
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;
    use crate::interval::Interval;
    type Iv = Interval<u32>;
    type Lapper<T> = crate::rust_lapper::Lapper<T, Interval<T>>;
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
        assert_eq!(None, lapper.super_find(15, 20).get(0));
        assert_eq!(None, lapper.seek(15, 20, &mut cursor).next());
    }

    // Test that a query start that hits an interval end returns no interval
    #[test]
    fn test_query_start_interval_stop() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(30, 35).next());
        assert_eq!(None, lapper.super_find(30, 35).get(0));
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
        assert_eq!(&expected, lapper.super_find(15, 25)[0]);
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
        assert_eq!(&expected, lapper.super_find(25, 35)[0]);
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
        assert_eq!(&expected, lapper.super_find(22, 27)[0]);
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
        assert_eq!(&expected, lapper.super_find(15, 35)[0]);
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
        assert_eq!(vec![&e1, &e2], lapper.super_find(8, 20));
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
        let found2 = lapper.super_find(8, 11);
        assert_eq!(found, vec![
            &Iv{start: 1, stop: 10, val: 0}, 
            &Iv{start: 9, stop: 11, val: 0},
            &Iv{start: 10, stop: 13, val: 0},
        ]);
        assert_eq!(found2, vec![
            &Iv{start: 1, stop: 10, val: 0}, 
            &Iv{start: 9, stop: 11, val: 0},
            &Iv{start: 10, stop: 13, val: 0},
        ]);
        let found = lapper.find(145, 151).collect::<Vec<&Iv>>();
        let found2 = lapper.super_find(145, 151);
        assert_eq!(found, vec![
            &Iv{start: 100, stop: 200, val: 0},
            &Iv{start: 111, stop: 160, val: 0},
            &Iv{start: 150, stop: 200, val: 0},
        ]);
        assert_eq!(found2, vec![
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
        let found2 = lapper.super_find(50, 55)[0];
        assert_eq!(found, Some(&e1));
        assert_eq!(found2, &e1);
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
        total = 0;
        for iv in lapper.iter() {
            total +=  lapper.super_find(iv.start, iv.stop).len();
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
        let found1 = lapper.super_find(28974798, 33141355);
        assert_eq!(found, vec![
            &Iv{start:28866309, stop: 33141404	, val: 0},
        ]);
        assert_eq!(found1, vec![
            &Iv{start:28866309, stop: 33141404	, val: 0},
        ]);
    }
}
