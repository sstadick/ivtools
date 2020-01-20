use crate::interval::Interval;

pub trait IvStore<'a, T>
where
    T: 'a,
{
    type IvIterator: Iterator<Item = &'a Interval<T>>;
    type IvSearchIterator: Iterator<Item = &'a Interval<T>>;
    type SelfU32: IvStore<'a, u32>;

    // functions that must be implemented
    fn new(intervals: Vec<Interval<T>>) -> Self;
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
