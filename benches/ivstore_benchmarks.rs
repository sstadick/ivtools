#[macro_use]
extern crate criterion;
extern crate ivtools;
extern crate rand;

use criterion::Criterion;
use ivtools::{Bits, Interval, IvStore, Lapper, ScAIList};
use rand::Rng;

type Iv = Interval<bool>;

fn randomi(imin: u32, imax: u32) -> u32 {
    let mut rng = rand::thread_rng();
    imin + rng.gen_range(0, imax - imin)
}

fn make_random(n: usize, range_max: u32, size_min: u32, size_max: u32) -> Vec<Iv> {
    let mut result = Vec::with_capacity(n);
    for _i in 0..n {
        let s = randomi(0, range_max);
        let e = s + randomi(size_min, size_max);
        result.push(Interval {
            start: s,
            stop: e,
            val: false,
        });
    }
    result
}

fn make_interval_set() -> (Vec<Iv>, Vec<Iv>) {
    //let n = 3_000_000;
    let n = 50_000;
    let chrom_size = 100_000_000;
    let min_interval_size = 500;
    let max_interval_size = 80000;
    let intervals = make_random(n, chrom_size, min_interval_size, max_interval_size);
    let other_intervals = make_random(n, 10 * chrom_size, 1, 2);
    (intervals, other_intervals)
}

pub fn query(c: &mut Criterion) {
    let s_size = 10;
    let (intervals, other_intervals) = make_interval_set();
    // Make Lapper intervals
    let mut bad_intervals: Vec<Iv> = intervals
        .iter()
        .map(|iv| Interval {
            start: iv.start,
            stop: iv.stop,
            val: true,
        })
        .collect();
    bad_intervals.push(Iv {
        start: 0,
        stop: 90_000_000,
        val: false,
    });

    let scailist = ScAIList::new(intervals.iter().map(|iv| Interval { ..*iv }).collect());
    let other_scailist = ScAIList::new(
        other_intervals
            .iter()
            .map(|iv| Interval { ..*iv })
            .collect(),
    );
    let bad_scailist = ScAIList::new(bad_intervals.iter().map(|iv| Interval { ..*iv }).collect());
    let bits = Bits::new(intervals.iter().map(|iv| Interval { ..*iv }).collect());
    // let other_bits = Bits::new(
    //     other_intervals
    //         .iter()
    //         .map(|iv| Interval { ..*iv })
    //         .collect(),
    // );
    let bad_bits = Bits::new(bad_intervals.iter().map(|iv| Interval { ..*iv }).collect());
    let lapper = Lapper::new(intervals);
    let other_lapper = Lapper::new(other_intervals);
    let bad_lapper = Lapper::new(bad_intervals);

    let mut comparison_group = c.benchmark_group("Bakeoff");
    comparison_group
        .sample_size(s_size)
        .bench_function("Lapper: find with 100% hit rate", |b| {
            b.iter(|| {
                for x in lapper.iter() {
                    lapper
                        .find(x.start, x.stop)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        });
    comparison_group
        .sample_size(s_size)
        .bench_function("Lapper: seek with 100% hit rate", |b| {
            b.iter(|| {
                let mut cursor = 0;
                for x in lapper.iter() {
                    lapper
                        .seek(x.start, x.stop, &mut cursor)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        });
    comparison_group.sample_size(s_size).bench_function(
        "Lapper: super_find with 100% hit rate",
        |b| {
            b.iter(|| {
                for x in lapper.iter() {
                    lapper.super_find(x.start, x.stop).len();
                }
            });
        },
    );
    comparison_group
        .sample_size(s_size)
        .bench_function("ScAIList: find with 100% hit rate", |b| {
            b.iter(|| {
                for x in scailist.iter() {
                    scailist
                        .find(x.start, x.stop)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        });
    comparison_group
        .sample_size(s_size)
        .bench_function("Bits: count with 100% hit rate", |b| {
            b.iter(|| {
                for x in lapper.iter() {
                    bits.count(x.start, x.stop);
                }
            });
        });

    comparison_group.sample_size(s_size).bench_function(
        "Lapper: find with below 100% hit rate",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    lapper
                        .find(x.start, x.stop)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        },
    );
    comparison_group.sample_size(s_size).bench_function(
        "Lapper: seek with below 100% hit rate",
        |b| {
            b.iter(|| {
                let mut cursor = 0;
                for x in other_lapper.iter() {
                    lapper
                        .seek(x.start, x.stop, &mut cursor)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        },
    );
    comparison_group.sample_size(s_size).bench_function(
        "Lapper: super find with below 100% hit rate",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    lapper.super_find(x.start, x.stop).len();
                }
            });
        },
    );
    comparison_group.sample_size(s_size).bench_function(
        "ScAIList: find with below 100% hit rate",
        |b| {
            b.iter(|| {
                for x in other_scailist.iter() {
                    scailist
                        .find(x.start, x.stop)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        },
    );
    comparison_group.sample_size(s_size).bench_function(
        "Bits: count with below 100% hit rate",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    bits.count(x.start, x.stop);
                }
            });
        },
    );

    comparison_group.sample_size(s_size).bench_function(
        "Lapper: find with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    bad_lapper
                        .find(x.start, x.stop)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        },
    );
    comparison_group.sample_size(s_size).bench_function(
        "Lapper: super find with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    bad_lapper.super_find(x.start, x.stop).len();
                }
            });
        },
    );
    comparison_group.sample_size(s_size).bench_function(
        "ScAIList: find with below 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in other_scailist.iter() {
                    bad_scailist
                        .find(x.start, x.stop)
                        .collect::<Vec<&Interval<bool>>>()
                        .len();
                }
            });
        },
    );
    comparison_group.sample_size(s_size).bench_function(
        "Bits: count with below 100% hit rate - chromosme spanning interval",
        |b| {
            b.iter(|| {
                for x in other_lapper.iter() {
                    bad_bits.count(x.start, x.stop);
                }
            });
        },
    );
    comparison_group.finish();
}
criterion_group!(benches, query);
criterion_main!(benches);
