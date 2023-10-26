use std::collections::HashMap;
use std::str::FromStr;

use pertable::Element;

use crate::atom::{Atom, PointChirality};
use crate::bond::{Bond, BondType};

enum AtomAttribute {
    Isotope,
    NImplicitHydrogens,
    Charge,
}

#[derive(Debug)]
pub enum MoleculeError {
    SmilesParseError(String),
}

#[derive(Debug, PartialEq, Eq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl FromStr for Molecule {
    type Err = MoleculeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut atoms = vec![];
        let mut bond = Bond::default();
        let mut bonds = vec![];
        let mut ring_closures = HashMap::new();

        let mut c_is_in_bracket = false;
        let mut atom_attribute = AtomAttribute::Isotope;
        let mut element_str = String::new();
        let mut root_atom = vec![];

        for c in s.chars() {
            if c_is_in_bracket {
                let atom: &mut Atom = atoms.last_mut().unwrap();
                if c == ']' {
                    c_is_in_bracket = false;
                    atom.element = match element_str.parse() {
                        Ok(element) => element,
                        Err(_) => {
                            return Err(MoleculeError::SmilesParseError(format!(
                                "{s} | invalid element {element_str}"
                            )))
                        }
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
                                "{s} | chirality error"
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
                || c == '.'
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
                    "{s} | invalid char {c}"
                )));
            }
        }

        Ok(Molecule { atoms, bonds })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::{Atom, PointChirality};
    use pretty_assertions::assert_eq;

    #[test]
    fn test_molecule_from_str() {
        let data = [
            (
                "C",
                Molecule {
                    atoms: vec![Atom {
                        element: Element::C,
                        isotope: None,
                        charge: 0,
                        delocalized: false,
                        n_implicit_hydrogens: None,
                        n_radical_electrons: None,
                        point_chirality: PointChirality::Undefined,
                    }],
                    bonds: vec![],
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
                        n_radical_electrons: None,
                        point_chirality: PointChirality::Undefined,
                    }],
                    bonds: vec![],
                },
            ),
            (
                "C-C",
                Molecule {
                    atoms: vec![
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
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![Bond {
                        i: 0,
                        j: 1,
                        bond_type: BondType::Single,
                    }],
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
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
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
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Default,
                        },
                    ],
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
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: Some(1),
                            n_radical_electrons: None,
                            point_chirality: PointChirality::CounterClockwise,
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
                        Atom {
                            element: Element::Br,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Default,
                        },
                    ],
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
                            element: Element::O,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::O,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Double,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Default,
                        },
                    ],
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
                            element: Element::F,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::F,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::F,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 4,
                            bond_type: BondType::Default,
                        },
                    ],
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
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                        Atom {
                            element: Element::H,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![Bond {
                        i: 0,
                        j: 1,
                        bond_type: BondType::Default,
                    }],
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
                            element: Element::C,
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
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 3,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 0,
                            j: 4,
                            bond_type: BondType::Default,
                        },
                    ],
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
                            element: Element::C,
                            isotope: None,
                            charge: 0,
                            delocalized: false,
                            n_implicit_hydrogens: None,
                            n_radical_electrons: None,
                            point_chirality: PointChirality::Undefined,
                        },
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 0,
                            j: 2,
                            bond_type: BondType::Default,
                        },
                    ],
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
                            element: Element::C,
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
                    ],
                    bonds: vec![
                        Bond {
                            i: 0,
                            j: 1,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 1,
                            j: 2,
                            bond_type: BondType::Default,
                        },
                        Bond {
                            i: 0,
                            j: 2,
                            bond_type: BondType::Double,
                        },
                        Bond {
                            i: 2,
                            j: 3,
                            bond_type: BondType::Default,
                        },
                    ],
                },
            ),
        ];

        for (smi, mol) in data {
            assert_eq!(Molecule::from_str(smi).unwrap(), mol);
        }
    }
}
