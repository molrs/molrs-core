use std::collections::HashMap;
use std::str::FromStr;

use pertable::Element;

use crate::atom::{Atom, PointChirality};
use crate::bond::{Bond, BondType};
use crate::utils::{deduplicate_vecs, get_index_of_duplicate};

enum AtomAttribute {
    Isotope,
    NImplicitHydrogens,
    Charge,
}

#[derive(Debug)]
pub enum MoleculeError {
    SmilesParseError(String),
    KekulizationError(String),
    BondOrderError(String),
    AssignImplicitHydrogensError(String),
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub rings: Option<Vec<Vec<usize>>>,
}

impl FromStr for Molecule {
    type Err = MoleculeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut mol = Molecule::parse_smi(s)?;

        match mol.perceive_implicit_hydrogens() {
            Ok(_) => (),
            Err(err) => {
                return Err(MoleculeError::AssignImplicitHydrogensError(format!(
                    "{s} | {:?}",
                    err
                )));
            }
        };
        mol.perceive_rings();
        mol = match mol.kekulized() {
            Ok(mol) => mol,
            Err(e) => {
                return Err(MoleculeError::KekulizationError(format!("{s} | {:?}", e)));
            }
        };
        mol.perceive_default_bonds();
        mol = mol.delocalized();
        // mol.perceive_stereo();  TODO

        Ok(mol)
    }
}

impl ToString for Molecule {
    fn to_string(&self) -> String {
        let atom = self.atoms[0];
        let mut atom_str = atom.to_string();
        if atom.element == Element::N && atom.delocalized && atom.n_implicit_hydrogens == Some(1) {
            atom_str = "[nH]".to_owned();
        };

        let mut smi = atom_str;
        let mut atom_str_indicies = vec![0];
        let mut ring_closure_count = 1;
        let mut ring_closure_bonds = vec![];

        for i in 1..self.atoms.len() {
            let atom = self.atoms[i];
            let mut atom_str = atom.to_string();

            if atom.element == Element::N
                && atom.delocalized
                && atom.n_implicit_hydrogens == Some(1)
            {
                atom_str = "[nH]".to_owned();
            };

            let mut neighbors = self.atom_neighbor_indicies(i);
            neighbors.retain(|neighbor| *neighbor < i);
            neighbors.sort();
            if neighbors.len() > 1 {
                for neighbor in neighbors.iter().take(neighbors.len() - 1) {
                    let ring_closure_bond = self.atoms_bond_between(i, *neighbor).unwrap();
                    ring_closure_bonds.push(ring_closure_bond);
                    let mut cursor = atom_str_indicies[ring_closure_bond.i + 1];
                    while ['(', ')'].contains(&smi.chars().nth(cursor).unwrap()) {
                        cursor -= 1;
                    }
                    let mut ring_closure_str = ring_closure_count.to_string();
                    ring_closure_count += 1;
                    smi.insert_str(cursor, &ring_closure_str);
                    if !(ring_closure_bond.bond_type == BondType::Single
                        || ring_closure_bond.bond_type == BondType::Delocalized
                        || ring_closure_bond.bond_type == BondType::Default)
                    {
                        ring_closure_str = String::from(char::from(ring_closure_bond.bond_type))
                            + &ring_closure_str;
                    }
                    atom_str += &ring_closure_str;
                    for index in atom_str_indicies[neighbor + 1..].iter_mut() {
                        *index += 1;
                    }
                }
            }

            neighbors = vec![*neighbors.last().unwrap()];

            if neighbors.contains(&(i - 1)) {
                atom_str_indicies.push(smi.len());
                smi += &atom_str;
            } else {
                let cursor = atom_str_indicies[neighbors[0]] + 1;
                let mut start_index = None;
                let mut parenthesis_depth = 0;
                for (i, c) in smi.chars().skip(cursor).enumerate() {
                    if c == '(' {
                        parenthesis_depth += 1;
                    } else if c == ')' {
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

                smi.insert(atom_str_indicies[start_index], '(');
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
                smi.insert(atom_str_indicies[bond.j], char::from(bond.bond_type));
                for index in atom_str_indicies[bond.j..].iter_mut() {
                    *index += 1;
                }
            }
        }

        smi
    }
}

impl Molecule {
    pub fn kekulized(&self) -> Result<Molecule, MoleculeError> {
        if self.rings.is_none() {
            return Err(MoleculeError::KekulizationError(format!(
                "{:?} | run perceive_rings before kekulization",
                self
            )));
        }

        let mut mol = self.clone();

        for bond in mol.bonds.iter_mut() {
            if bond.bond_type == BondType::Delocalized {
                bond.bond_type = BondType::Default;
            }
        }

        let mut conjugated_rings = vec![];
        for ring in mol.rings.as_ref().unwrap() {
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
                    let atom = mol.atoms[index];
                    let maximum_valence = atom.element.valence(atom.charge).unwrap();
                    let bond_order =
                        mol.atom_explicit_valence(index) + atom.n_implicit_hydrogens.unwrap();
                    if mol.atom_n_double_bonds(index) == 0 && bond_order < maximum_valence {
                        for path in contiguous_paths.iter_mut() {
                            if self
                                .atoms_bond_between(index, *path.last().unwrap())
                                .is_some()
                            {
                                path.push(index);
                            } else if self
                                .atoms_bond_between(index, *path.first().unwrap())
                                .is_some()
                            {
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
                        return Err(MoleculeError::KekulizationError(format!("{:?}", self)));
                    } else if path.len() % 2 == 1 {
                        needs_kekulization.push(true);
                    } else {
                        needs_kekulization.push(false);
                        for i in 0..(path.len() / 2) {
                            mol.atoms_bond_between_mut(path[i * 2], path[i * 2 + 1])
                                .unwrap()
                                .bond_type = BondType::Double;
                            mol.atoms[i * 2].delocalized = false;
                            mol.atoms[i * 2 + 1].delocalized = false;
                        }
                    }
                }
            }
            if needs_kekulization.is_empty() {
                break;
            }
            if needs_kekulization.iter().all(|val| *val) {
                return Err(MoleculeError::KekulizationError(format!("{:?}", self)));
            }
        }

        for bond in mol.bonds.iter_mut() {
            if bond.bond_type == BondType::Default {
                bond.bond_type = BondType::Single;
            }
        }

        for atom in mol.atoms.iter_mut() {
            atom.delocalized = false;
        }

        Ok(mol)
    }

    pub fn delocalized(&self) -> Molecule {
        let mut mol = self.clone();

        for ring in self.rings.as_ref().unwrap().iter().rev() {
            if ring
                .iter()
                .all(|index| self.atom_needs_delocalization(*index))
            {
                for index in ring {
                    mol.atoms[*index].delocalized = true;
                }
                mol.atoms_bond_between_mut(*ring.first().unwrap(), *ring.last().unwrap())
                    .unwrap()
                    .bond_type = BondType::Delocalized;
                for window in ring.windows(2) {
                    mol.atoms_bond_between_mut(window[0], window[1])
                        .unwrap()
                        .bond_type = BondType::Delocalized;
                }
            }
        }

        mol
    }

    pub fn atom_bonds(&self, index: usize) -> Vec<&Bond> {
        self.bonds
            .iter()
            .filter(|bond| bond.i == index || bond.j == index)
            .collect()
    }

    pub fn atoms_bond_between(&self, i: usize, j: usize) -> Option<&Bond> {
        self.bonds
            .iter()
            .find(|&bond| (bond.i == i && bond.j == j) || (bond.i == j && bond.j == i))
    }

    pub fn atoms_bond_between_mut(&mut self, i: usize, j: usize) -> Option<&mut Bond> {
        self.bonds
            .iter_mut()
            .find(|bond| (bond.i == i && bond.j == j) || (bond.i == j && bond.j == i))
    }

    pub fn atom_neighbor_indicies(&self, index: usize) -> Vec<usize> {
        let mut neighbor_indicies = vec![];
        for bond in &self.bonds {
            if bond.i == index {
                neighbor_indicies.push(bond.j);
            } else if bond.j == index {
                neighbor_indicies.push(bond.i);
            }
        }

        neighbor_indicies
    }

    pub fn atom_explicit_valence(&self, index: usize) -> u8 {
        self.atom_bonds(index)
            .iter()
            .map(|bond| f64::from(bond.bond_type))
            .sum::<f64>() as u8
    }

    fn parse_smi(smi: &str) -> Result<Molecule, MoleculeError> {
        let mut atoms = vec![];
        let mut bond = Bond::default();
        let mut bonds = vec![];
        let mut ring_closures = HashMap::new();

        let mut c_is_in_bracket = false;
        let mut atom_attribute = AtomAttribute::Isotope;
        let mut element_str = String::new();
        let mut root_atom = vec![];

        for c in smi.chars() {
            if c_is_in_bracket {
                let atom: &mut Atom = atoms.last_mut().unwrap();
                if c == ']' {
                    c_is_in_bracket = false;
                    atom.element = match element_str.parse() {
                        Ok(element) => element,
                        Err(_) => {
                            return Err(MoleculeError::SmilesParseError(format!(
                                "{smi} | invalid element {element_str}"
                            )))
                        }
                    };
                    if element_str.chars().nth(0).unwrap().is_lowercase() {
                        atom.delocalized = true;
                    };
                    element_str = String::new();
                } else if c == 'H' {
                    if element_str.is_empty() {
                        element_str.push('H');
                    } else {
                        atom.n_implicit_hydrogens = Some(1);
                        atom_attribute = AtomAttribute::NImplicitHydrogens;
                    }
                } else if c == '+' {
                    atom.charge = 1;
                    atom_attribute = AtomAttribute::Charge;
                } else if c == '-' {
                    atom.charge = -1;
                    atom_attribute = AtomAttribute::Charge;
                } else if c == '@' {
                    match atom.point_chirality {
                        PointChirality::Undefined => {
                            atom.point_chirality = PointChirality::CounterClockwise;
                        }
                        PointChirality::CounterClockwise => {
                            atom.point_chirality = PointChirality::Clockwise;
                        }
                        _ => {
                            return Err(MoleculeError::SmilesParseError(format!(
                                "{smi} | chirality error"
                            )));
                        }
                    }
                } else if c.is_alphabetic() {
                    element_str.push(c);
                } else if c.is_numeric() {
                    let c_as_digit = c.to_digit(10).unwrap();
                    match atom_attribute {
                        AtomAttribute::Isotope => match atom.isotope {
                            None => {
                                atom.isotope = Some(c_as_digit as u16);
                            }
                            Some(isotope) => atom.isotope = Some(isotope * 10 + c_as_digit as u16),
                        },
                        AtomAttribute::Charge => {
                            atom.charge *= c_as_digit as i8;
                        }
                        AtomAttribute::NImplicitHydrogens => {
                            atom.n_implicit_hydrogens = Some(c_as_digit as u8);
                        }
                    };
                }
            } else if c == '[' {
                c_is_in_bracket = true;
                atom_attribute = AtomAttribute::Isotope;
                atoms.push(Atom::default());
                if !root_atom.is_empty() {
                    bond.i = root_atom.pop().unwrap();
                    bond.j = atoms.len() - 1;
                    bonds.push(bond);
                    bond = Bond::default();
                }
                root_atom.push(atoms.len() - 1);
            } else if c == '-'
                || c == '/'
                || c == '\\'
                || c == ':'
                || c == '='
                || c == '#'
                || c == '$'
            {
                bond.bond_type = BondType::try_from(c).unwrap();
            } else if c == '%' {
                unimplemented!()
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap() as usize;
                if let std::collections::hash_map::Entry::Vacant(e) =
                    ring_closures.entry(c_as_digit)
                {
                    e.insert(atoms.len() - 1);
                } else {
                    let ring_closure_bond = Bond {
                        i: *ring_closures.get(&c_as_digit).unwrap(),
                        j: atoms.len() - 1,
                        bond_type: bond.bond_type,
                    };
                    bond.bond_type = BondType::Default;
                    bonds.push(ring_closure_bond);
                    ring_closures.remove(&c_as_digit);
                }
            } else if c == '(' {
                root_atom.push(*root_atom.last().unwrap());
            } else if c == ')' {
                root_atom.pop();
            } else if c == 'b'
                || c == 'c'
                || c == 'n'
                || c == 'o'
                || c == 'p'
                || c == 's'
                || c == 'B'
                || c == 'C'
                || c == 'N'
                || c == 'O'
                || c == 'P'
                || c == 'S'
                || c == 'F'
                || c == 'I'
                || c == '*'
            {
                atoms.push(Atom {
                    element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                    isotope: None,
                    charge: 0,
                    delocalized: c.is_lowercase(),
                    n_implicit_hydrogens: None,
                    n_radical_electrons: None,
                    point_chirality: crate::atom::PointChirality::Undefined,
                });
                if !root_atom.is_empty() {
                    bond.i = root_atom.pop().unwrap();
                    bond.j = atoms.len() - 1;
                    bonds.push(bond);
                    bond = Bond::default();
                }
                root_atom.push(atoms.len() - 1);
            } else if c == 'l' && atoms.last().unwrap().element == Element::C {
                atoms.last_mut().unwrap().element = Element::Cl;
            } else if c == 'r' && atoms.last().unwrap().element == Element::B {
                atoms.last_mut().unwrap().element = Element::Br;
            } else {
                return Err(MoleculeError::SmilesParseError(format!(
                    "{smi} | invalid char {c}"
                )));
            }
        }

        Ok(Molecule {
            atoms,
            bonds,
            rings: None,
        })
    }

    fn perceive_implicit_hydrogens(&mut self) -> Result<(), MoleculeError> {
        let bond_orders: Vec<u8> = (0..self.atoms.len())
            .map(|index| self.atom_explicit_valence(index))
            .collect();
        let vec_n_double_bonds: Vec<u8> = (0..self.atoms.len())
            .map(|index| self.atom_n_double_bonds(index))
            .collect();

        for (i, atom) in self.atoms.iter_mut().enumerate() {
            let mut maximum_valence = atom.element.valence(atom.charge).unwrap();
            let n_double_bonds = vec_n_double_bonds[i];
            match atom.element {
                Element::P => {
                    if n_double_bonds == 0 {
                    } else if n_double_bonds == 1 {
                        maximum_valence += 2;
                    } else {
                        return Err(MoleculeError::BondOrderError(
                            "P has more than 1 double bond".to_owned(),
                        ));
                    }
                }
                Element::S => {
                    if n_double_bonds == 0 {
                    } else if n_double_bonds < 3 {
                        maximum_valence += 2 * n_double_bonds;
                    } else {
                        return Err(MoleculeError::BondOrderError(
                            "S has more than 2 double bonds".to_owned(),
                        ));
                    }
                }
                Element::Cl => {
                    if n_double_bonds == 0 {
                    } else if n_double_bonds < 4 {
                        maximum_valence += 2 * n_double_bonds;
                    } else {
                        return Err(MoleculeError::BondOrderError(
                            "Cl has more than 3 double bonds".to_owned(),
                        ));
                    }
                }
                _ => (),
            };
            let bond_order = bond_orders[i];
            if bond_order > maximum_valence {
                return Err(MoleculeError::BondOrderError(
                    "explicit valence is higher than maximum allowed valence".to_owned(),
                ));
            }
            let mut n_implicit_hydrogens = maximum_valence - bond_order;
            if atom.n_implicit_hydrogens.is_none() {
                if atom.delocalized && n_implicit_hydrogens > 0 {
                    n_implicit_hydrogens -= 1;
                }
                atom.n_implicit_hydrogens = Some(n_implicit_hydrogens);
                atom.n_radical_electrons = Some(0);
            } else {
                atom.n_radical_electrons =
                    Some(n_implicit_hydrogens - atom.n_implicit_hydrogens.unwrap());
            }
        }

        Ok(())
    }

    fn perceive_rings(&mut self) {
        let mut paths = vec![];
        let mut closed_paths = vec![];
        for neighbor_index in self.atom_neighbor_indicies(0) {
            paths.push(vec![0, neighbor_index]);
        }

        let mut counter = 0;
        while counter < 20 && !paths.is_empty() {
            counter += 1;

            let mut new_paths = vec![];
            for path in paths.iter_mut() {
                if path.is_empty() {
                    continue;
                }
                let last_atom_index = *path.last().unwrap();
                let second_to_last_atom_index = path.iter().rev().nth(1).unwrap();
                let mut neighbor_indices = self.atom_neighbor_indicies(last_atom_index);
                neighbor_indices
                    .retain(|neighbor_index| neighbor_index != second_to_last_atom_index);
                if neighbor_indices.is_empty() {
                    *path = vec![];
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
                    if let Some(index_of_duplicate) = get_index_of_duplicate(path) {
                        closed_paths.push(path[(index_of_duplicate)..path.len() - 1].to_owned());
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
        closed_paths = deduplicate_vecs(closed_paths);
        closed_paths.sort_by_key(|closed_loop| -(closed_loop.len() as isize));

        self.rings = Some(closed_paths);
    }

    fn perceive_default_bonds(&mut self) {
        for bond in self.bonds.iter_mut() {
            if bond.bond_type != BondType::Default {
                continue;
            }
            if self.atoms.get(bond.i).unwrap().delocalized
                && self.atoms.get(bond.j).unwrap().delocalized
            {
                bond.bond_type = BondType::Delocalized;
            } else {
                bond.bond_type = BondType::Single;
            }
        }
    }

    fn atom_n_double_bonds(&self, index: usize) -> u8 {
        self.atom_bonds(index)
            .iter()
            .filter(|bond| bond.bond_type == BondType::Double)
            .count() as u8
    }

    fn atom_needs_kekulization(&self, index: usize) -> bool {
        let atom = self.atoms[index];
        atom.delocalized || self.atom_n_double_bonds(index) > 0
    }

    fn atom_needs_delocalization(&self, index: usize) -> bool {
        let atom = self.atoms[index];
        self.atom_n_double_bonds(index) > 0
            || (self.atom_bonds(index).len() as u8 + atom.n_implicit_hydrogens.unwrap()) as i8
                + atom.charge
                < 4
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::{Atom, PointChirality};
    use pretty_assertions::assert_eq;

    #[test]
    fn playground() {
        let mol = Molecule::from_str("[nH]1cncc1").unwrap();
        dbg!(&mol);
    }

    #[test]
    fn test_molecule_from_and_to_str() {
        let data = [
            (
                "C",
                Molecule {
                    atoms: vec![Atom {
                        element: Element::C,
                        isotope: None,
                        charge: 0,
                        delocalized: false,
                        n_implicit_hydrogens: Some(4),
                        n_radical_electrons: Some(0),
                        point_chirality: PointChirality::Undefined,
                    }],
                    bonds: vec![],
                    rings: Some(vec![]),
                },
            ),
            (
                "[13CH3+]",
                Molecule {
                    atoms: vec![Atom {
                        element: Element::C,
                        isotope: Some(13),
                        charge: 1,
                        delocalized: false,
                        n_implicit_hydrogens: Some(3),
                        n_radical_electrons: Some(0),
                        point_chirality: PointChirality::Undefined,
                    }],
                    bonds: vec![],
                    rings: Some(vec![]),
                },
            ),
            (
                "CC",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![Bond {
                        i: 0,
                        j: 1,
                        bond_type: BondType::Single,
                    }],
                    rings: Some(vec![]),
                },
            ),
            (
                "C#N",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::N,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![Bond {
                        i: 0,
                        j: 1,
                        bond_type: BondType::Triple,
                    }],
                    rings: Some(vec![]),
                },
            ),
            (
                "C[CH]C",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(1),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Single,
                        },
                    ],
                    rings: Some(vec![]),
                },
            ),
            (
                "F[C@H](Cl)Br",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::F,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::CounterClockwise,
                        },
                        Atom {
                            element: Element::Cl,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::Br,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Single,
                        },
                    ],
                    rings: Some(vec![]),
                },
            ),
            (
                "CC(=O)O",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::O,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::O,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Double,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Single,
                        },
                    ],
                    rings: Some(vec![]),
                },
            ),
            (
                "CC(F)(F)F",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::F,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::F,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::F,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 4,
                            bond_type: BondType::Single,
                        },
                    ],
                    rings: Some(vec![]),
                },
            ),
            (
                "C[H]",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::H,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![Bond {
                        i: 0,
                        j: 1,
                        bond_type: BondType::Single,
                    }],
                    rings: Some(vec![]),
                },
            ),
            (
                "C(C(C)C)C",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(2),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 0,
                            j: 4,
                            bond_type: BondType::Single,
                        },
                    ],
                    rings: Some(vec![]),
                },
            ),
            (
                "C1CC1",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(2),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(2),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(2),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 0,
                            j: 2,
                            bond_type: BondType::Single,
                        },
                    ],
                    rings: Some(vec![vec![0, 1, 2]]),
                },
            ),
            (
                "C1CC=1C",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(2),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(0),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(3),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Single,
                        },
                        Bond {
                            i: 0,
                            j: 2,
                            bond_type: BondType::Double,
                        },
                        Bond {
                            i: 2,
                            j: 3,
                            bond_type: BondType::Single,
                        },
                    ],
                    rings: Some(vec![vec![0, 1, 2]]),
                },
            ),
            (
                "c1ccccc1",
                Molecule {
                    atoms: vec![
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: true,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: true,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: true,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: true,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: true,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: true,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: Some(0),
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Delocalized,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Delocalized,
                        },
                        Bond {
                            i: 2,
                            j: 3,
                            bond_type: BondType::Delocalized,
                        },
                        Bond {
                            i: 3,
                            j: 4,
                            bond_type: BondType::Delocalized,
                        },
                        Bond {
                            i: 4,
                            j: 5,
                            bond_type: BondType::Delocalized,
                        },
                        Bond {
                            i: 0,
                            j: 5,
                            bond_type: BondType::Delocalized,
                        },
                    ],
                    rings: Some(vec![vec![0, 1, 2, 3, 4, 5]]),
                },
            ),
        ];

        for (smi, mol) in data {
            dbg!(&smi);
            let mol_from_smi = Molecule::from_str(smi).unwrap();
            assert_eq!(mol_from_smi, mol);
            assert_eq!(mol_from_smi.to_string(), smi);
        }
    }

    #[test]
    fn test_kekulize() {
        let mol = Molecule::from_str("c1ccccn1").unwrap();
        assert_eq!("C1=CC=CC=N1", mol.kekulized().unwrap().to_string());

        let mol = Molecule::from_str("[nH]1cccc1").unwrap();
        assert_eq!("N1C=CC=C1", mol.kekulized().unwrap().to_string());

        let mol = Molecule::from_str("c1ccc2ccccn12").unwrap();
        assert_eq!("C1=CC=C2C=CC=CN12", mol.kekulized().unwrap().to_string());
    }

    #[test]
    fn test_delocalized() {
        let mol = Molecule::from_str("C1C=CC=CC=1").unwrap();
        assert_eq!("c1ccccc1", mol.to_string());

        let mol = Molecule::from_str("N1C=CN=C1").unwrap();
        assert_eq!("[nH]1ccnc1", mol.to_string());

        let mol = Molecule::from_str("C12=CC=CN1N=CS2").unwrap();
        assert_eq!("c12cccn1ncs2", mol.to_string());
    }

    #[test]
    fn test_atom_bonds() {
        let mol = Molecule::from_str("CC(C)C").unwrap();
        assert_eq!(
            mol.atom_bonds(1),
            vec![
                &Bond {
                    i: 0,
                    j: 1,
                    bond_type: BondType::Single
                },
                &Bond {
                    i: 1,
                    j: 2,
                    bond_type: BondType::Single
                },
                &Bond {
                    i: 1,
                    j: 3,
                    bond_type: BondType::Single
                }
            ]
        );
    }

    #[test]
    fn test_atoms_bond_between() {
        let mol = Molecule::from_str("CC#N").unwrap();
        assert_eq!(
            mol.atoms_bond_between(1, 2),
            Some(&Bond {
                i: 1,
                j: 2,
                bond_type: BondType::Triple
            })
        );
    }

    #[test]
    fn test_atom_neighbor_indicies() {
        let mol = Molecule::from_str("CCC").unwrap();
        assert_eq!(mol.atom_neighbor_indicies(1), vec![0, 2]);
    }

    #[test]
    fn test_atom_explicit_valence() {
        let mol = Molecule::from_str("CS(=O)C").unwrap();
        assert_eq!(mol.atom_explicit_valence(1), 4);
    }
}
