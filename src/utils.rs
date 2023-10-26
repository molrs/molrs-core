use std::collections::HashSet;

pub fn get_index_of_duplicate(v: &[usize]) -> Option<usize> {
    let mut hashset = HashSet::new();
    for item in v {
        if hashset.contains(item) {
            return Some(v.iter().position(|val| val == item).unwrap());
        } else {
            hashset.insert(*item);
        }
    }

    None
}

pub fn deduplicate_vecs(mut closed_loops: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut sorted_clone = closed_loops.clone();
    sorted_clone
        .iter_mut()
        .for_each(|closed_loop| closed_loop.sort());
    let mut loops_to_drop = HashSet::new();
    for i in (1..sorted_clone.len()).rev() {
        let closed_loop = sorted_clone.get(i).unwrap();
        for other_loop in &sorted_clone[0..i] {
            if closed_loop == other_loop {
                loops_to_drop.insert(i);
            }
        }
    }
    let mut loops_to_drop: Vec<usize> = loops_to_drop.into_iter().collect();
    loops_to_drop.sort();
    for i in loops_to_drop.iter().rev() {
        closed_loops.remove(*i);
    }

    closed_loops
}
