use std::collections::HashSet;

use petgraph::{Graph, Undirected};

use crate::{
    atom::Atom,
    bond::BondType,
    utils::{deduplicate_nested_vec, get_duplicate_element, traverse_paths},
};

#[derive(Debug, Clone, Default)]
pub struct Molecule {
    pub graph: Graph<Atom, BondType, Undirected, usize>,
    pub rings: Vec<Vec<usize>>,
}

impl Molecule {
    pub fn perceive_rings(&mut self) {
        if self.graph.node_count() < 3 {
            return;
        }

        let mut paths = vec![];
        for neighbor in self.graph.neighbors(0.into()) {
            paths.push(vec![0, neighbor.index()]);
        }
        fn condition(path: &mut Vec<usize>) -> Result<Vec<usize>, ()> {
            if let Some(duplicate_element) = get_duplicate_element(path) {
                for (i, index) in path.iter().enumerate() {
                    if *index == duplicate_element {
                        return Ok(path[(i + 1)..].to_owned());
                    }
                }
            }

            Err(())
        }
        let mut paths = deduplicate_nested_vec(&traverse_paths(&self.graph, &paths, condition));
        paths.sort_by_key(|path| path.len());
        self.rings = paths;
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
        let mut num_rad_electrons = vec![];
        for (i, atom) in self.graph.node_weights().enumerate() {
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
            if atom.needs_update {
                num_imp_hs.push(num_imp_h);
            } else {
                num_rad_electrons.push(num_imp_h - atom.num_imp_h)
            }
        }
        num_imp_hs.reverse();
        num_rad_electrons.reverse();
        for atom in self.graph.node_weights_mut() {
            if atom.needs_update {
                atom.needs_update = false;
                atom.num_imp_h = num_imp_hs.pop().unwrap();
            } else {
                atom.num_rad_electron = num_rad_electrons.pop().unwrap();
            }
        }

        Ok(())
    }

    pub fn kekulized(&self) -> Result<Molecule, String> {
        let mut mol = self.clone();

        let mut conjugated_rings = vec![];
        for ring in &self.rings {
            if ring.iter().all(|index| {
                self.graph.node_weight((*index).into()).unwrap().aromatic
                    || self
                        .graph
                        .edges((*index).into())
                        .any(|edge| *edge.weight() == BondType::Double)
            }) {
                conjugated_rings.push(ring.clone())
            }
        }

        let mut needs_kekulization = vec![true];
        while needs_kekulization.iter().any(|val| *val) {
            needs_kekulization = vec![];
            for conjugated_ring in &conjugated_rings {
                let mut contiguous_paths: Vec<Box<Vec<usize>>> = vec![];
                for index in conjugated_ring.iter() {
                    let atom = mol.graph.node_weight((*index).into()).unwrap();
                    let maximum_valence = atom.atomic_symbol.maximum_valence(atom.charge, true);
                    let bond_order = mol
                        .graph
                        .edges((*index).into())
                        .map(|x| x.weight().to_float())
                        .sum::<f64>() as u8
                        + atom.num_imp_h;
                    if !mol
                        .graph
                        .edges((*index).into())
                        .any(|edge| *edge.weight() == BondType::Double)
                        && bond_order <= maximum_valence
                    {
                        for path in contiguous_paths.iter_mut() {
                            if self
                                .graph
                                .find_edge((*index).into(), (*path.last().unwrap()).into())
                                .is_some()
                            {
                                path.push(*index);
                            } else if self
                                .graph
                                .find_edge((*index).into(), (*path.first().unwrap()).into())
                                .is_some()
                            {
                                let mut new_path = vec![*index];
                                new_path.extend(path.iter());
                                *path = Box::new(new_path);
                            }
                        }
                        if !contiguous_paths.iter().any(|path| path.contains(index)) {
                            contiguous_paths.push(Box::new(vec![*index]));
                        }
                    }
                }

                for path in contiguous_paths.iter() {
                    if path.len() == 1 {
                        return Err("kekulization failed - isolated aromatic atom".to_owned());
                    } else if path.len() % 2 == 1 {
                        needs_kekulization.push(true);
                    } else {
                        needs_kekulization.push(false);
                        for i in 0..(path.len() / 2) {
                            mol.graph.update_edge(
                                path[i * 2].into(),
                                path[i * 2 + 1].into(),
                                BondType::Double,
                            );
                        }
                    }
                }
            }
            if needs_kekulization.is_empty() {
                break;
            }
            if needs_kekulization.iter().all(|val| *val) {
                return Err("kekulization failed - all paths have odd length".to_owned());
            }
        }

        for edge in mol.graph.edge_weights_mut() {
            if *edge == BondType::Aromatic {
                *edge = BondType::Single;
            }
        }

        for atom in mol.graph.node_weights_mut() {
            atom.aromatic = false;
        }

        Ok(mol)
    }

    // pub fn aromatized(&self) -> Result<Molecule, String> {
    //     let mol = self.clone();

    //     Ok(mol)
    // }

    fn unique_rings(&self) -> Vec<Vec<usize>> {
        let mut unique_rings = vec![];
        let mut unique_atoms = HashSet::new();
        for ring in &self.rings {
            if !ring.iter().all(|index| unique_atoms.contains(index)) {
                ring.iter().for_each(|index| {
                    unique_atoms.insert(*index);
                });
                unique_rings.push(ring.clone());
            }
        }
        unique_rings.reverse();

        unique_rings
    }

    pub fn coordinates_2d(&self) -> Vec<Option<[f64; 2]>> {
        let mut unique_rings = self.unique_rings();
        let mut _unique_chains = self.unique_chains();

        let mut coords: Vec<Option<[f64; 2]>> =
            (0..self.graph.node_count()).map(|_| None).collect();
        coords[0] = Some([0.0, 0.0]);
        let n_neighbors = self.graph.neighbors(0.into()).count();
        if n_neighbors == 1 {
            coords[self.graph.neighbors(0.into()).next().unwrap().index()] = Some([0.0, 1.0]);
        } else if n_neighbors == 2 {
            let mut neighbors_iter = self.graph.neighbors(0.into());
            coords[neighbors_iter.next().unwrap().index()] = Some([0.0, 1.0]);
            coords[neighbors_iter.next().unwrap().index()] = Some([-0.866, -0.5]);
        } else if n_neighbors == 3 {
            unimplemented!();
        } else {
            unimplemented!();
        }

        for (i, coord) in coords.iter_mut().enumerate().skip(1) {
            if coord.is_some() {
                continue;
            }

            let mut rings_to_remove = vec![];
            for (j, ring) in unique_rings.iter().enumerate() {
                if ring.contains(&i) {
                    rings_to_remove.push(j);
                    unimplemented!();
                }
            }
            for j in rings_to_remove.iter().rev() {
                unique_rings.remove(*j);
            }
            if coord.is_some() {
                continue;
            }

            let n_neighbors = self.graph.neighbors(i.into()).count();
            if n_neighbors == 1 {}
        }

        // coords.iter().map(|coord| coord.unwrap()).collect()
        coords
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_unique_chains() {
        // let smi = "C";
        // let smi = "CCC";
        // let smi = "C(C)C";
        let smi = "CC(C)C";
        // let smi = "CCCCC";
        // let smi = "c12ncccc1[nH]cc2";
        // let smi = "C1C(C)CCN1C";
        // let smi = "C1(CC2)CC2CCC1";
        let mol = crate::from::smiles(smi).unwrap();
        // dbg!(&mol);
        dbg!(mol.unique_chains());
    }

    #[test]
    fn test_coordinates_2d() {
        // let smi = "CCC";
        let smi = "C(C)C";
        // let smi = "c12ncccc1[nH]cc2";
        // let smi = "C1(CC2)CC2CCC1";
        let mol = crate::from::smiles(smi).unwrap();
        // dbg!(&mol);
        dbg!(mol.coordinates_2d());
    }
}
