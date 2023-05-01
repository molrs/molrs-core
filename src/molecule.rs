use petgraph::{Graph, Undirected};

use crate::{atom::Atom, bond::BondType, utils::get_duplicate_element};

#[derive(Debug, Clone, Default)]
pub struct Molecule {
    pub graph: Graph<Atom, BondType, Undirected, usize>,
    pub rings: Vec<Box<Vec<usize>>>,
}

impl Molecule {
    pub fn perceive_rings(&mut self) {
        if self.graph.node_count() < 3 {
            return;
        }

        let mut paths = vec![];
        for neighbor in self.graph.neighbors(0.into()) {
            paths.push(Box::new(vec![0, neighbor.index()]));
        }

        while !paths.is_empty() {
            let mut new_paths = vec![];
            for path in paths.iter_mut() {
                if path.len() == 0 {
                    continue;
                }
                let last_atom_index = path.last().unwrap();
                let second_to_last_atom_index = path.iter().rev().nth(1).unwrap();
                let neighbors: Vec<usize> = self
                    .graph
                    .neighbors((*last_atom_index).into())
                    .map(|neighbor| neighbor.index())
                    .filter(|neighbor_index| neighbor_index != second_to_last_atom_index)
                    .collect();
                let neighbors_len = neighbors.len();
                if neighbors_len == 0 {
                    *path = Box::<std::vec::Vec<usize>>::default();
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
            for i in (0..paths.len()).rev() {
                if paths[i].len() == 0 {
                    paths.remove(i);
                }
            }
            for new_path in new_paths {
                paths.push(new_path);
            }

            for path in paths.iter_mut() {
                if let Some(duplicate_element) = get_duplicate_element(path) {
                    for (i, index) in path.iter().enumerate() {
                        if *index == duplicate_element {
                            self.rings.push(Box::new(path[(i + 1)..].to_owned()));
                            break;
                        }
                    }
                    *path = Box::<std::vec::Vec<usize>>::default();
                }
            }
        }

        let mut rings = self.rings.clone();
        rings.iter_mut().for_each(|ring| ring.sort());
        let mut drop_indices = vec![];
        for i in (1..rings.len()).rev() {
            let ring = rings.get(i).unwrap();
            for other_ring in &rings[0..i] {
                if ring == other_ring {
                    drop_indices.push(i);
                }
            }
        }

        for i in (0..self.rings.len()).rev() {
            if drop_indices.contains(&i) {
                self.rings.remove(i);
            }
        }
    }

    pub fn perceive_default_bonds(&mut self) {
        let is_aromatic_vec: Vec<bool> = (0..self.graph.edge_count())
            .map(|i| self.graph.edge_endpoints(i.into()).unwrap())
            .map(|tup| {
                self.graph.node_weight(tup.0).unwrap().aromatic
                    && self.graph.node_weight(tup.1).unwrap().aromatic
            })
            .collect();
        for (bond_type, is_aromatic) in self.graph.edge_weights_mut().zip(is_aromatic_vec) {
            if bond_type == &BondType::Default {
                if is_aromatic {
                    *bond_type = BondType::Aromatic;
                } else {
                    *bond_type = BondType::Single;
                }
            }
        }
    }

    pub fn perceive_implicit_hydrogens(&mut self) -> Result<(), String> {
        let mut num_imp_hs = vec![];
        for (i, atom) in self.graph.node_weights().enumerate() {
            if atom.needs_update {
                let maximum_valence = atom.atomic_symbol.maximum_valence(atom.charge, true);
                let bond_order = self
                    .graph
                    .edges(i.into())
                    .map(|x| x.weight().to_float())
                    .sum::<f64>() as u8;
                if bond_order > maximum_valence {
                    return Err("bond order is too high".to_owned());
                }
                let mut num_imp_h = maximum_valence - bond_order;
                while num_imp_h > atom.atomic_symbol.maximum_valence(atom.charge, false) {
                    num_imp_h -= 2;
                }
                num_imp_hs.push(num_imp_h);
            }
        }
        num_imp_hs.reverse();
        for atom in self.graph.node_weights_mut() {
            if atom.needs_update {
                atom.needs_update = false;
                atom.num_imp_h = num_imp_hs.pop().unwrap();
            }
        }

        Ok(())
    }

    // pub fn kekulized(&self) -> Result<Molecule, String> {
    //     let mut mol = self.clone();

    //     Ok(mol)
    // }

    // pub fn aromatized(&self) -> Result<Molecule, String> {
    //     let mut mol = self.clone();

    //     Ok(mol)
    // }
}
