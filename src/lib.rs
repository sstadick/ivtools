pub mod interval;
pub mod ivstore;
pub mod rust_lapper;
pub mod scailist;

// Test the methods against eachother
#[cfg(test)]
mod tests {
    use crate::interval::Interval;
    use crate::ivstore::IvStore;
    use crate::rust_lapper::Lapper;
    use crate::scailist::ScAIList;

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

    #[test]
    fn test_non_overlapping() {
        let data1 = setup_nonoverlapping();
        let data2 = setup_nonoverlapping();
        let lapper = Lapper::new(data1);
        let scailist = ScAIList::new(data2);

        let mut result1: Vec<&Iv> = lapper
            .iter()
            .flat_map(|iv| lapper.find(iv.start, iv.stop).collect::<Vec<&Iv>>())
            .collect();
        let mut result2: Vec<&Iv> = scailist
            .iter()
            .flat_map(|iv| scailist.find(iv.start, iv.stop).collect::<Vec<&Iv>>())
            .collect();
        assert_eq!(result1.sort(), result2.sort())
    }

    #[test]
    fn test_overlapping() {
        let data1 = setup_overlapping();
        let data2 = setup_overlapping();
        let lapper = Lapper::new(data1);
        let scailist = ScAIList::new(data2);

        let mut result1: Vec<&Iv> = lapper
            .iter()
            .flat_map(|iv| lapper.find(iv.start, iv.stop).collect::<Vec<&Iv>>())
            .collect();
        let mut result2: Vec<&Iv> = scailist
            .iter()
            .flat_map(|iv| scailist.find(iv.start, iv.stop).collect::<Vec<&Iv>>())
            .collect();
        assert_eq!(result1.sort(), result2.sort())
    }
    #[test]
    fn test_bad() {
        let data1 = setup_bad();
        let data2 = setup_bad();
        let lapper = Lapper::new(data1);
        let scailist = ScAIList::new(data2);

        let mut result1: Vec<&Iv> = lapper
            .iter()
            .flat_map(|iv| lapper.find(iv.start, iv.stop).collect::<Vec<&Iv>>())
            .collect();
        let mut result2: Vec<&Iv> = scailist
            .iter()
            .flat_map(|iv| scailist.find(iv.start, iv.stop).collect::<Vec<&Iv>>())
            .collect();
        assert_eq!(result1.sort(), result2.sort())
    }
}
