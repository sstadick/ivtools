pub trait IvStore<'a, T, I>
where
    T: 'a,
    I: IntervalLike<T>,
{
    type IvIterator: Iterator<Item = &'a I>;
    type IvSearchIterator: Iterator<Item = &'a I>;
    type SelfU32: IvStore<'a, u32, I>;

    // functions that must be implemented
    fn new(intervals: Vec<I>) -> Self;
    fn len(&'a self) -> usize;
    fn is_empty(&'a self) -> bool;
    fn iter(&'a self) -> Self::IvIterator;
    fn merge_overlaps(self) -> Self::SelfU32;
    fn find(&'a self, start: u32, stop: u32) -> Self::IvSearchIterator;
    // fn depth(&'a self) -> Self::IvIterator;
    // fn coverage(&'a self) -> u32;
    // fn intersect(&'a self, other: &Self) -> u32;
    // fn union(&'a self, other: &Self) -> u32;
}

pub trait IntervalLike<T>: Ord + Eq + PartialEq + PartialOrd {
    fn new(start: u32, stop: u32, val: T) -> Self;
    fn intersect(&self, other: &Self) -> u32 {
        std::cmp::min(self.stop(), other.stop())
            .checked_sub(std::cmp::max(self.start(), other.start()))
            .unwrap_or(0)
    }
    fn overlap(&self, start: u32, stop: u32) -> bool {
        self.start() < stop && self.stop() > start
    }
    fn start(&self) -> u32;
    fn stop(&self) -> u32;
    fn val(&self) -> &T;
    fn set_stop(&mut self, new: u32);
    fn set_start(&mut self, new: u32);
    fn set_val(&mut self, new: T);
}
