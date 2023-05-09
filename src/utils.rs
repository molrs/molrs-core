use std::collections::HashSet;

use crate::{atom::Atom, bond::BondType};
use petgraph::{Graph, Undirected};

pub fn get_duplicate_element(vector: &Vec<usize>) -> Option<usize> {
    let mut set = HashSet::new();
    for value in vector {
        if set.contains(value) {
            return Some(*value);
        } else {
            set.insert(value);
        }
    }

    None
}

pub fn traverse_paths(
    graph: &Graph<Atom, BondType, Undirected, usize>,
    initial_paths: &[Vec<usize>],
    condition: impl Fn(&mut Vec<usize>) -> Result<Vec<usize>, ()>,
) -> Vec<Vec<usize>> {
    let mut paths = initial_paths.to_owned();
    let mut paths_to_return = vec![];

    while !paths.is_empty() {
        let mut new_paths = vec![];
        for path in paths.iter_mut() {
            if path.is_empty() {
                continue;
            }
            let last_atom_index = path.last().unwrap();
            let second_to_last_atom_index = path.iter().rev().nth(1).unwrap();
            let neighbors: Vec<usize> = graph
                .neighbors((*last_atom_index).into())
                .map(|neighbor| neighbor.index())
                .filter(|neighbor_index| neighbor_index != second_to_last_atom_index)
                .collect();
            let neighbors_len = neighbors.len();
            if neighbors_len == 0 {
                *path = vec![];
            } else if neighbors_len == 1 {
                path.push(neighbors[0]);
            } else {
                for neighbor in &neighbors[1..] {
                    let mut new_path = path.clone();
                    new_path.push(*neighbor);
                    new_paths.push(new_path);
                }
                path.push(neighbors[0]);
            }
        }
        for new_path in new_paths {
            paths.push(new_path);
        }
        for path in paths.iter_mut() {
            if !path.is_empty() {
                if let Ok(path_to_return) = condition(path) {
                    paths_to_return.push(path_to_return);
                    *path = vec![];
                }
            }
        }
        for i in (0..paths.len()).rev() {
            if paths[i].is_empty() {
                paths.remove(i);
            }
        }
    }

    paths_to_return
}

pub fn deduplicate_nested_vec(nested_vec: &[Vec<usize>], mode: &str) -> Vec<Vec<usize>> {
    let mut nested_vec = nested_vec.to_owned();
    let mut tmp_nested_vec = nested_vec.clone();
    if mode == "sort" {
        tmp_nested_vec.iter_mut().for_each(|vec| vec.sort());
    } else if mode.is_empty() {
    } else {
        unimplemented!("not a valid deduplication mode");
    }
    let mut drop_indices = vec![];
    for i in (1..tmp_nested_vec.len()).rev() {
        let ring = tmp_nested_vec.get(i).unwrap();
        for other_ring in &tmp_nested_vec[0..i] {
            if ring == other_ring {
                drop_indices.push(i);
            }
        }
    }
    for i in (0..tmp_nested_vec.len()).rev() {
        if drop_indices.contains(&i) {
            nested_vec.remove(i);
        }
    }

    nested_vec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_duplicate_elements() {
        assert_eq!(get_duplicate_element(&vec![0, 1, 2, 3, 4, 0]).unwrap(), 0);
        assert_eq!(get_duplicate_element(&vec![0, 1, 2, 3, 4]), Option::None);
    }
}
