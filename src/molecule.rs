use std::str::FromStr;

use pertable::Element;

use crate::atom::Atom;
use crate::bond::{Bond, BondType};

#[derive(Debug)]
pub enum MoleculeError {
    SmilesParseError(String),
}

#[derive(Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl FromStr for Molecule {
    type Err = MoleculeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut atoms = vec![];
        let mut bonds = vec![];

        let mut atom = Atom::default();
        let mut bond = Bond::default();

        let mut c_is_in_bracket = false;

        for c in s.chars() {
            if c_is_in_bracket {
                // modify atom based on what c is
                if c == ']' {
                    c_is_in_bracket = false;
                }
            } else if c == '[' {
                c_is_in_bracket = true;
                atoms.push(atom);
                atom = Atom::default();
                bonds.push(bond);
                bond = Bond::default();
            } else if c == '-'
                || c == '/'
                || c == '\\'
                || c == ':'
                || c == '='
                || c == '#'
                || c == '$'
                || c == '.'
            {
                bond.bond_type = BondType::try_from(c).unwrap();
            } else if c == '%' || c.is_numeric() {
                // handle ring closures
                unimplemented!()
            } else if c == '(' || c == ')' {
                // handle parentheses
                unimplemented!()
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
                atom.element = Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap();
                atom.delocalized = c.is_lowercase();
                atoms.push(atom);
                atom = Atom::default();
                bonds.push(bond);
                bond = Bond::default();
            } else if c == 'l' && atoms.last().unwrap().element == Element::C {
                atoms.last_mut().unwrap().element = Element::Cl;
            } else if c == 'r' && atoms.last().unwrap().element == Element::B {
                atoms.last_mut().unwrap().element = Element::Br;
            } else {
                return Err(MoleculeError::SmilesParseError(format!(
                    "{s} | invalid char {c}"
                )));
            }
        }

        bonds.pop();

        Ok(Molecule { atoms, bonds })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::{Atom, PointChirality};

    #[test]
    fn test_single_atom_molecule_from_str() {
        // add some bracket atoms later
        let smiles = ["*", "C", "Cl", "A", "Zn"];
        let expected_fails = [false, false, false, true, true];
        let atoms = [
            Atom {
                element: Element::Any,
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: PointChirality::Undefined,
            },
            Atom {
                element: Element::C,
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: PointChirality::Undefined,
            },
            Atom {
                element: Element::Cl,
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: PointChirality::Undefined,
            },
        ];

        for ((smi, expected_fail), atom) in smiles.into_iter().zip(expected_fails).zip(atoms) {
            let mol = Molecule::from_str(smi);
            if expected_fail {
                assert!(mol.is_err());
            } else {
                assert!(mol.is_ok());
                assert_eq!(*mol.unwrap().atoms.first().unwrap(), atom);
            }
        }
    }
}
