use std::{cmp::Ordering, str::FromStr};

use crate::{
    atom::Atom,
    bond::{Bond, BondType},
    smiles::{SmilesParseError, SmilesParser},
    utils::{deduplicate_closed_loops, get_duplicate},
};

pub struct KekulizationError;

pub struct BondOrderError;

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    rings: Vec<Vec<usize>>,
}

impl FromStr for Molecule {
    type Err = SmilesParseError;

    /// Parses a SMILES to a Molecule.
    ///
    /// First SMILES are sliced into its component atom_strs, bond_strs, and
    /// ring_closure_strs. Then atom_strs are parsed to atoms and bond_strs +
    /// ring_closure_strs are parsed to bonds.
    ///
    /// The molecule is then processed by finding rings, converting Default bond
    /// types to Single or Delocalized bond types, kekulizing Delocalized bonds,
    /// and perceiving implicit hydrogens and radical electrons.
    fn from_str(smi: &str) -> Result<Self, Self::Err> {
        let smiles_parser: SmilesParser = smi.parse()?;
        let atoms = smiles_parser.atoms()?;
        let mut bonds = smiles_parser.bonds()?;
        bonds.sort();

        let mut mol = Molecule {
            atoms,
            bonds,
            rings: vec![],
        };
        mol.perceive_rings();
        mol.perceive_default_bonds();
        mol = match mol.kekulized() {
            Ok(mol) => mol,
            Err(_) => {
                return Err(SmilesParseError {
                    details: format!("{} | kekulization failed", smi),
                })
            }
        };
        let is_bracket_atom: Vec<bool> = smiles_parser
            .atom_strs
            .iter()
            .map(|atom_str| atom_str.starts_with('['))
            .collect();
        match mol.perceive_implicit_hydrogens(is_bracket_atom) {
            Ok(_) => (),
            Err(_) => {
                return Err(SmilesParseError {
                    details: format!("{} | bond order error", smi),
                })
            }
        }
        mol = mol.delocalized();
        // perceive stereochemistry

        Ok(mol)
    }
}

impl ToString for Molecule {
    fn to_string(&self) -> String {
        let mut smi = self.atoms.first().unwrap().to_string();
        let mut atom_str_indicies = vec![0];
        let mut ring_closure_count = 1;
        let mut ring_closure_bonds = vec![];

        for i in 1..self.atoms.len() {
            let mut atom_str = self.atoms.get(i).unwrap().to_string();

            let mut neighbors = self.neighbor_indices(i);
            neighbors.retain(|neighbor| *neighbor < i);
            neighbors.sort();
            if neighbors.len() > 1 {
                for j in 0..neighbors.len() - 1 {
                    let neighbor = neighbors[j];

                    let ring_closure_bond = self.bond_between(neighbor, i).unwrap();
                    ring_closure_bonds.push(ring_closure_bond);
                    let mut cursor = *atom_str_indicies.get(ring_closure_bond.atom_i + 1).unwrap();
                    while ['(', ')'].contains(&smi.chars().nth(cursor).unwrap()) {
                        cursor -= 1;
                    }
                    let mut ring_closure_str = match ring_closure_count.cmp(&10) {
                        Ordering::Greater => format!("%{}", ring_closure_count),
                        Ordering::Less => format!("{}", ring_closure_count),
                        Ordering::Equal => format!("%{}", ring_closure_count),
                    };
                    ring_closure_count += 1;
                    smi.insert_str(cursor, &ring_closure_str);
                    if !(ring_closure_bond.bond_type == BondType::Single
                        || ring_closure_bond.bond_type == BondType::Delocalized
                        || ring_closure_bond.bond_type == BondType::Default)
                    {
                        ring_closure_str =
                            String::from(ring_closure_bond.bond_type.to_char()) + &ring_closure_str;
                    }
                    atom_str += &ring_closure_str;
                    for index in atom_str_indicies[neighbor + 1..].iter_mut() {
                        *index += 1;
                    }

                    neighbors.remove(j);
                }
            }

            if neighbors.contains(&(i - 1)) {
                atom_str_indicies.push(smi.len());
                smi += &atom_str;
            } else {
                let cursor = atom_str_indicies.get(neighbors[0]).unwrap() + 1;
                let mut start_index = None;
                let mut parenthesis_depth = 0;
                for (i, c) in smi.chars().skip(cursor).enumerate() {
                    if c == '(' {
                        parenthesis_depth += 1;
                    }
                    if c == ')' {
                        parenthesis_depth -= 1;
                    }
                    if parenthesis_depth > 0 {
                        continue;
                    }
                    if let Some(index) = atom_str_indicies
                        .iter()
                        .position(|index| *index == cursor + i)
                    {
                        start_index = Some(index);
                        break;
                    }
                }

                let start_index = start_index.unwrap();

                smi.insert(*atom_str_indicies.get(start_index).unwrap(), '(');
                for index in atom_str_indicies[start_index..].iter_mut() {
                    *index += 1;
                }

                smi.push(')');
                atom_str_indicies.push(smi.len());
                smi += &atom_str;
            }
        }

        for bond in &self.bonds {
            if !(bond.bond_type == BondType::Default
                || bond.bond_type == BondType::Single
                || bond.bond_type == BondType::Delocalized
                || ring_closure_bonds.contains(&bond))
            {
                smi.insert(
                    *atom_str_indicies.get(bond.atom_j).unwrap(),
                    bond.bond_type.to_char(),
                );
                for index in atom_str_indicies[bond.atom_j..].iter_mut() {
                    *index += 1;
                }
            }
        }

        smi
    }
}

impl Molecule {
    /// Returns the explicit valence for an atom. Does not include implicit
    /// hydrogens.
    pub fn explicit_valence(&self, index: usize) -> u8 {
        self.bonds_to_atom(index)
            .iter()
            .map(|bond| bond.bond_type.to_float())
            .sum::<f64>() as u8
    }

    /// Returns the atom indices of the neighbors of an atom.
    pub fn neighbor_indices(&self, index: usize) -> Vec<usize> {
        let mut neighbor_indices = vec![];
        for bond in &self.bonds {
            if bond.atom_i == index {
                neighbor_indices.push(bond.atom_j);
            } else if bond.atom_j == index {
                neighbor_indices.push(bond.atom_i);
            }
        }
        neighbor_indices.sort();

        neighbor_indices
    }

    /// Returns references to all bonds to an atom.
    pub fn bonds_to_atom(&self, index: usize) -> Vec<&Bond> {
        self.bonds
            .iter()
            .filter(|bond| bond.atom_i == index || bond.atom_j == index)
            .collect()
    }

    /// Returns a reference to the bond between atom_i and atom_j.
    pub fn bond_between(&self, atom_i: usize, atom_j: usize) -> Option<&Bond> {
        self.bonds.iter().find(|&bond| {
            (bond.atom_i == atom_i && bond.atom_j == atom_j)
                || (bond.atom_i == atom_j && bond.atom_j == atom_i)
        })
    }

    /// Returns a mutable reference to the bond between atom_i and atom_j.
    pub fn bond_between_mut(&mut self, atom_i: usize, atom_j: usize) -> Option<&mut Bond> {
        self.bonds.iter_mut().find(|bond| {
            (bond.atom_i == atom_i && bond.atom_j == atom_j)
                || (bond.atom_i == atom_j && bond.atom_j == atom_i)
        })
    }

    /// Returns a kekulized clone of the Molecule.
    pub fn kekulized(&self) -> Result<Molecule, KekulizationError> {
        let mut mol = self.clone();

        let mut conjugated_rings = vec![];
        for ring in &self.rings {
            if ring
                .iter()
                .all(|index| self.atom_needs_kekulization(*index))
            {
                conjugated_rings.push(ring.clone());
            }
        }

        let mut needs_kekulization = vec![true];
        while needs_kekulization.iter().any(|val| *val) {
            needs_kekulization = vec![];
            for conjugated_ring in &conjugated_rings {
                let mut contiguous_paths: Vec<Vec<usize>> = vec![];
                for index in conjugated_ring {
                    let index = *index;
                    let atom = mol.atoms.get(index).unwrap();
                    let maximum_valence = atom.element.maximum_valence(atom.charge, true);
                    let bond_order = mol.explicit_valence(index) + atom.num_implicit_hydrogens;
                    if !mol.atom_has_double_bond(index) && bond_order <= maximum_valence {
                        for path in contiguous_paths.iter_mut() {
                            if self.bond_between(index, *path.last().unwrap()).is_some() {
                                path.push(index);
                            } else if self.bond_between(index, *path.first().unwrap()).is_some() {
                                let mut new_path = vec![index];
                                new_path.extend(path.iter());
                                *path = new_path;
                            }
                        }
                        if !contiguous_paths.iter().any(|path| path.contains(&index)) {
                            contiguous_paths.push(vec![index]);
                        }
                    }
                }

                for path in contiguous_paths.iter() {
                    if path.len() == 1 {
                        return Err(KekulizationError);
                    } else if path.len() % 2 == 1 {
                        needs_kekulization.push(true);
                    } else {
                        needs_kekulization.push(false);
                        for i in 0..(path.len() / 2) {
                            mol.bond_between_mut(path[i * 2], path[i * 2 + 1])
                                .unwrap()
                                .bond_type = BondType::Double;
                        }
                    }
                }
            }
            if needs_kekulization.is_empty() {
                break;
            }
            if needs_kekulization.iter().all(|val| *val) {
                return Err(KekulizationError);
            }
        }

        for bond in mol.bonds.iter_mut() {
            if bond.bond_type == BondType::Delocalized {
                bond.bond_type = BondType::Single;
            }
        }

        for atom in mol.atoms.iter_mut() {
            atom.delocalized = false;
        }

        Ok(mol)
    }

    /// Returns a delocalized clone of the Molecule.
    pub fn delocalized(&self) -> Molecule {
        let mut mol = self.clone();

        for ring in self.rings.iter().rev() {
            if ring
                .iter()
                .all(|index| self.atom_needs_delocalization(*index))
            {
                dbg!(&ring);
                for index in ring {
                    mol.atoms[*index].delocalized = true;
                }
                mol.bond_between_mut(*ring.first().unwrap(), *ring.last().unwrap())
                    .unwrap()
                    .bond_type = BondType::Delocalized;
                for window in ring.windows(2) {
                    mol.bond_between_mut(window[0], window[1])
                        .unwrap()
                        .bond_type = BondType::Delocalized;
                }
            }
        }

        mol
    }

    /// Traverses the Molecule to find all rings.
    fn perceive_rings(&mut self) {
        let mut paths = vec![];
        let mut closed_loops = vec![];
        for neighbor_index in self.neighbor_indices(0) {
            paths.push(vec![0, neighbor_index]);
        }

        while !paths.is_empty() {
            let mut new_paths = vec![];
            for path in paths.iter_mut() {
                if path.is_empty() {
                    continue;
                }
                let last_atom_index = *path.last().unwrap();
                let second_to_last_atom_index = path.iter().rev().nth(1).unwrap();
                let mut neighbor_indices = self.neighbor_indices(last_atom_index);
                neighbor_indices
                    .retain(|neighbor_index| neighbor_index != second_to_last_atom_index);
                if neighbor_indices.is_empty() {
                    *path = vec![];
                } else if neighbor_indices.len() == 1 {
                    path.push(neighbor_indices[0]);
                } else {
                    for neighbor_index in &neighbor_indices[1..] {
                        let mut new_path = path.clone();
                        new_path.push(*neighbor_index);
                        new_paths.push(new_path);
                    }
                    path.push(neighbor_indices[0]);
                }
            }
            for new_path in new_paths {
                paths.push(new_path);
            }
            for path in paths.iter_mut() {
                if !path.is_empty() {
                    if let Some(duplicate) = get_duplicate(path) {
                        let index = path.iter().position(|value| *value == duplicate).unwrap();
                        closed_loops.push(path[(index + 1)..].to_owned());
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
        closed_loops = deduplicate_closed_loops(closed_loops);
        closed_loops.sort_by_key(|closed_loop| -(closed_loop.len() as isize));

        self.rings = closed_loops;
    }

    /// Converts Default bond types to Single or Delocalized bond types.
    fn perceive_default_bonds(&mut self) {
        for bond in self.bonds.iter_mut() {
            if bond.bond_type != BondType::Default {
                continue;
            }
            if self.atoms.get(bond.atom_i).unwrap().delocalized
                && self.atoms.get(bond.atom_j).unwrap().delocalized
            {
                bond.bond_type = BondType::Delocalized;
            } else {
                bond.bond_type = BondType::Single;
            }
        }
    }

    /// Computes the number of implicit hydrogens when implicit hydrogens are
    /// not specified or the number of radical electrons when implicit hydrogens
    /// are specified.
    fn perceive_implicit_hydrogens(
        &mut self,
        is_bracket_atom: Vec<bool>,
    ) -> Result<(), BondOrderError> {
        let bond_orders: Vec<u8> = (0..self.atoms.len())
            .map(|index| self.explicit_valence(index))
            .collect();
        for (i, atom) in self.atoms.iter_mut().enumerate() {
            let maximum_valence = atom.element.maximum_valence(atom.charge, true);
            let bond_order = bond_orders[i];
            if bond_order > maximum_valence {
                return Err(BondOrderError);
            }
            let mut num_imp_hs = maximum_valence - bond_order;
            while num_imp_hs > atom.element.maximum_valence(atom.charge, false) {
                num_imp_hs -= 2;
            }
            if !is_bracket_atom[i] {
                atom.num_implicit_hydrogens = num_imp_hs;
            } else {
                atom.num_radical_electrons = num_imp_hs - atom.num_implicit_hydrogens;
            }
        }

        Ok(())
    }

    /// Returns whether an atom participates in a double bond or not.
    fn atom_has_double_bond(&self, index: usize) -> bool {
        return self
            .bonds_to_atom(index)
            .iter()
            .any(|bond| bond.bond_type == BondType::Double);
    }

    /// Returns whether an atom is delocalized or it participates in a double
    /// bond.
    fn atom_needs_kekulization(&self, index: usize) -> bool {
        let atom = self.atoms.get(index).unwrap();
        atom.delocalized || self.atom_has_double_bond(index)
    }

    fn atom_needs_delocalization(&self, index: usize) -> bool {
        let atom = self.atoms.get(index).unwrap();
        self.atom_has_double_bond(index)
            || (self.explicit_valence(index) + atom.num_implicit_hydrogens) as i8 + atom.charge < 4
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_string_from_molecule() {
        let smiles = [
            "CC",
            "CCCCCC",
            "C[18O-]",
            "C[CH]C",
            "C=C",
            "CC(C)C",
            "C(=O)C",
            "CS(=O)(=O)C",
            "C(C(C)C)C",
            "C1CC1",
            "C1CCCC1",
            "c12ccccc1[nH]cc2",
            "c2ccc1ccccc1c2",
            "c13[nH]ncc1[C@H]2C[C@H]2C3",
            "NC(Cn1c3c(c(C(F)(F)F)n1)[C@H]2C[C@H]2C3(F)F)=O",
            "CC(C)(C#Cc1nc(c(cc1)c2c3c(c(cc2)Cl)c(nn3CC(F)(F)F)[NH]S(=O)(=O)C)[C@H](Cc4cc(cc(c4)F)F)[NH]C(=O)Cn7c6c([C@H]5C[C@H]5C6(F)F)c(n7)C(F)(F)F)S(=O)(=O)C",
        ];
        for smi in smiles {
            let mol: Molecule = smi.parse().unwrap();
            assert_eq!(mol.to_string(), smi);
        }
    }

    #[test]
    fn test_molecule_from_str() {
        // let smi = "C1CC1";
        // let smi = "CC2CC2C";
        // let smi = "C12=CC=CC=C1NC=C2";
        // let smi = "C1=CC=C2C=CC=CC2=C1";
        // let smi = "c1ccc1";
        let smi = "C1=CC=C1";
        // let smi = "Cc1cc(C)c1";
        // let smi = "CC(C)(C)C";
        let mol: Molecule = smi.parse().unwrap();

        dbg!(&mol);
    }
}
