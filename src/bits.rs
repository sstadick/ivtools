use crate::interval::Interval;

pub struct Bits {
    /// Sorted list of start positions
    starts: Vec<u32>,
    /// Sorted list of end positions,
    stops: Vec<u32>,
}

impl Bits {
    pub fn new<T>(intervals: Vec<Interval<T>>) -> Self {
        let (mut starts, mut stops): (Vec<_>, Vec<_>) =
            intervals.iter().map(|x| (x.start, x.stop)).unzip();
        starts.sort();
        stops.sort();
        Bits {
            starts: starts,
            stops: stops,
        }
    }
    #[inline]
    pub fn bsearch_seq(key: u32, elems: &[u32]) -> usize {
        if elems[0] > key {
            return 0;
        }
        let mut high = elems.len();
        let mut low = 0;

        while high - low > 1 {
            let mid = (high + low) / 2;
            if elems[mid] < key {
                low = mid;
            } else {
                high = mid;
            }
        }
        high
    }

    /// Count all intervals that overlap start .. stop. This performs two binary search in order to
    /// find all the excluded elements, and then deduces the intersection from there. See
    /// [BITS](https://arxiv.org/pdf/1208.3407.pdf) for more details.
    /// ```
    /// use ivtools::{Bits, Interval};
    /// let bits = Bits::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<bool>>>());
    /// assert_eq!(bits.count(5, 11), 2);
    /// ```
    #[inline]
    pub fn count(&self, start: u32, stop: u32) -> usize {
        let len = self.starts.len();
        let mut first = Self::bsearch_seq(start, &self.stops);
        let last = Self::bsearch_seq(stop, &self.starts);
        //println!("{}/{}", start, stop);
        //println!("pre start found in stops: {}: {}", first, self.stops[first]);
        //println!("pre stop found in starts: {}", last);
        //while last < len && self.starts[last] == stop {
        //last += 1;
        //}
        while first < len && self.stops[first] == start {
            first += 1;
        }
        let num_cant_after = len - last;
        let result = len - first - num_cant_after;
        //println!("{:#?}", self.starts);
        //println!("{:#?}", self.stops);
        //println!("start found in stops: {}", first);
        //println!("stop found in starts: {}", last);
        result
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod test {
    use super::*;
    use crate::interval::Interval;
    use crate::rust_lapper::Lapper;
    use crate::ivstore::IvStore;

    
        type Iv = Interval<u32>;
        fn setup_nonoverlapping() -> Vec<Iv> {
            (0..100_000)
                .step_by(20)
                .map(|x| Iv {
                    start: x,
                    stop: x + 10,
                    val: 0,
                })
                .collect()
        }
        fn setup_overlapping() -> Vec<Iv> {
            (0..100_000)
                .step_by(10)
                .map(|x| Iv {
                    start: x,
                    stop: x + 15,
                    val: 0,
                })
                .collect()
        }
        fn setup_bad() -> Vec<Iv> {
            vec![
                Iv {
                    start: 70,
                    stop: 120,
                    val: 0,
                }, // max_len = 50
                Iv {
                    start: 10,
                    stop: 15,
                    val: 0,
                },
                Iv {
                    start: 10,
                    stop: 15,
                    val: 0,
                }, // exact overlap
                Iv {
                    start: 12,
                    stop: 15,
                    val: 0,
                }, // inner overlap
                Iv {
                    start: 14,
                    stop: 16,
                    val: 0,
                }, // overlap end
                Iv {
                    start: 40,
                    stop: 45,
                    val: 0,
                },
                Iv {
                    start: 50,
                    stop: 55,
                    val: 0,
                },
                Iv {
                    start: 60,
                    stop: 65,
                    val: 0,
                },
                Iv {
                    start: 68,
                    stop: 71,
                    val: 0,
                }, // overlap start
                Iv {
                    start: 70,
                    stop: 75,
                    val: 0,
                },
            ]
        }
    
        // Test that a query stop that hits an interval start returns no interval
        #[test]
        fn test_query_stop_interval_start() {
            let data1 = setup_nonoverlapping();
            let data2 = setup_nonoverlapping();
            let lapper = Lapper::new(data1);
            let bits = Bits::new(data2);
            assert_eq!(lapper.find(15, 20).count(), bits.count(15, 20));
        }
    
        // Test that a query start that hits an interval end returns no interval
        #[test]
        fn test_query_start_interval_stop() {
            let data1 = setup_nonoverlapping();
            let data2 = setup_nonoverlapping();
            let lapper = Lapper::new(data1);
            let bits = Bits::new(data2);

            assert_eq!(lapper.find(30, 35).count(), bits.count(30, 35));
        }
    
        // Test that a query that overlaps the start of an interval returns that interval
        #[test]
        fn test_query_overlaps_interval_start() {
            let data1 = setup_nonoverlapping();
            let data2 = setup_nonoverlapping();
            let lapper = Lapper::new(data1);
            let bits = Bits::new(data2);

            assert_eq!(lapper.find(15, 25).count(), bits.count(15, 25));
        }
    
        // Test that a query that overlaps the stop of an interval returns that interval
        #[test]
        fn test_query_overlaps_interval_stop() {
            let data1 = setup_nonoverlapping();
            let data2 = setup_nonoverlapping();
            let lapper = Lapper::new(data1);
            let bits = Bits::new(data2);
            assert_eq!(lapper.find(25, 35).count(), bits.count(25, 35));
        }
    
        // Test that a query that is enveloped by interval returns interval
        #[test]
        fn test_interval_envelops_query() {
            let data1 = setup_nonoverlapping();
            let data2 = setup_nonoverlapping();
            let lapper = Lapper::new(data1);
            let bits = Bits::new(data2);
            assert_eq!(lapper.find(22, 27).count(), bits.count(22, 27));
        }
    
        // Test that a query that envolops an interval returns that interval
        #[test]
        fn test_query_envolops_interval() {
            let data1 = setup_nonoverlapping();
            let data2 = setup_nonoverlapping();
            let lapper = Lapper::new(data1);
            let bits = Bits::new(data2);
            assert_eq!(lapper.find(15, 35).count(), bits.count(15, 35));
        }
    
        #[test]
        fn test_overlapping_intervals() {
            let data2 = setup_overlapping();
            let bits = Bits::new(data2);
            assert_eq!(bits.count(8, 20), 2);
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
            let bits = Bits::new(data1);
            assert_eq!(bits.count(8, 11), 3);
            assert_eq!(bits.count(145, 151), 3);
        }
        // BUG TESTS - these are tests that came from real life
    

    
        // Test that if lower_bound puts us before the first match, we still return a match
        #[test]
        fn test_find_over_behind_first_match() {
            let data1 = setup_bad();
            let data2 = setup_bad();
            let lapper = Lapper::new(data1);
            let bits = Bits::new(data2);
            assert_eq!(lapper.find(50, 55).count(), bits.count(50,55));
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
            let bits = Bits::new(data);
            assert_eq!(bits.count(28974798, 33141355), 1);
        }
}
