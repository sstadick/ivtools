pub struct Bits {
    /// Sorted list of start positions
    starts: Vec<u32>,
    /// Sorted list of end positions,
    stops: Vec<u32>,
}

impl Bits {
    pub fn new(starts: Vec<u32>, stops: Vec<u32>) -> Self {
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
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, stop: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<bool>>>());
    /// assert_eq!(lapper.count(5, 11), 2);
    /// ```
    #[inline]
    pub fn count(&self, start: u32, stop: u32) -> usize {
        let len = self.intervals.len();
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
