use std::str::FromStr;

use crate::periodic_table::AtomicSymbol;

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Chirality {
    #[default]
    Undefined,
    Clockwise,
    CounterClockwise,
}

pub struct ChiralityParseError;

impl FromStr for Chirality {
    type Err = ChiralityParseError;

    fn from_str(chirality: &str) -> Result<Self, Self::Err> {
        if chirality.is_empty() {
            Ok(Chirality::Undefined)
        } else if chirality == "@@" {
            Ok(Chirality::Clockwise)
        } else if chirality == "@" {
            Ok(Chirality::CounterClockwise)
        } else {
            return Result::Err(ChiralityParseError);
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Atom {
    pub index: usize,
    pub atomic_symbol: AtomicSymbol,
    pub atomic_number: u8,
    pub isotope: u16,
    pub charge: i8,
    pub aromatic: bool,
    pub num_imp_h: u8,
    pub chirality: Chirality,
    pub num_rad_electron: usize,
    pub ring_sizes: Vec<usize>,
    pub needs_update: bool,
}

enum NextValue {
    Isotope,
    Charge(bool),
    NumImpH,
}

impl Atom {
    pub fn from_str(atom_str: &str, index: usize) -> Result<Atom, String> {
        if atom_str.starts_with('[') {
            let mut atomic_symbol = String::new();
            let mut atom = Atom::default();
            let mut next_value = NextValue::Isotope;
            for c in atom_str.chars() {
                if c == '[' {
                } else if c == ']' {
                    break;
                } else if c.is_numeric() {
                    let c_as_int = c.to_digit(10).unwrap();
                    match next_value {
                        NextValue::Isotope => {
                            atom.isotope = atom.isotope * 10 + c_as_int as u16;
                        }
                        NextValue::Charge(is_positive) => {
                            if is_positive {
                                atom.charge = c_as_int as i8;
                            } else {
                                atom.charge = -(c_as_int as i8);
                            }
                        }
                        NextValue::NumImpH => {
                            atom.num_imp_h = c_as_int as u8;
                        }
                    }
                } else if c.is_alphabetic() {
                    let len_atomic_symbol = atomic_symbol.chars().count();
                    if len_atomic_symbol == 0 {
                        atomic_symbol.push(c);
                        atom.aromatic = c.is_lowercase();
                    } else if len_atomic_symbol == 1 && c != 'H' {
                        atomic_symbol.push(c);
                    } else if c == 'H' {
                        atom.num_imp_h = 1;
                        next_value = NextValue::NumImpH;
                    } else {
                    }
                } else if c == '@' {
                    atom.chirality = match atom.chirality {
                        Chirality::Undefined => Chirality::CounterClockwise,
                        Chirality::CounterClockwise => Chirality::Clockwise,
                        _ => return Err(format!("chirality error in atom_str {}", &atom_str)),
                    }
                } else if c == '+' {
                    atom.charge = 1;
                    next_value = NextValue::Charge(true);
                } else if c == '-' {
                    atom.charge = -1;
                    next_value = NextValue::Charge(false);
                } else {
                }
            }
            let atomic_symbol = match AtomicSymbol::from_str(&atomic_symbol) {
                Ok(atomic_symbol) => atomic_symbol,
                Err(err) => return Err(err),
            };
            atom.index = index;
            atom.atomic_number = atomic_symbol.atomic_number();
            atom.atomic_symbol = atomic_symbol;
            Ok(atom)
        } else {
            let mut atomic_symbol = String::new();
            let mut atom = Atom::default();
            for c in atom_str.chars() {
                if !c.is_alphabetic() && c != '*' {
                    break;
                }
                if atomic_symbol.chars().count() == 0 {
                    atom.aromatic = c.is_lowercase();
                }
                atomic_symbol.push(c);
            }
            let atomic_symbol = match AtomicSymbol::from_str(&atomic_symbol) {
                Ok(atomic_symbol) => atomic_symbol,
                Err(err) => return Err(err),
            };
            atom.index = index;
            atom.atomic_number = atomic_symbol.atomic_number();
            atom.atomic_symbol = atomic_symbol;
            atom.needs_update = true;
            Ok(atom)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::periodic_table::AtomicSymbol;

    use super::{Atom, Chirality};

    #[test]
    fn test_atom_from_str() {
        let atom_strs = [
            "*",
            "C",
            "Cl",
            "c",
            "[H]",
            "[13C]",
            "[H+]",
            "[O-2]",
            "[ClH]",
            "[C@H]",
            "[2H+]",
            "[13CH3]",
            "[13CH3+]",
            "[13C@@H+]",
        ];
        let atoms = [
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::Star,
                atomic_number: 0,
                isotope: 0,
                charge: 0,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: true,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::C,
                atomic_number: 6,
                isotope: 0,
                charge: 0,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: true,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::Cl,
                atomic_number: 17,
                isotope: 0,
                charge: 0,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: true,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::C,
                atomic_number: 6,
                isotope: 0,
                charge: 0,
                aromatic: true,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: true,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::H,
                atomic_number: 1,
                isotope: 0,
                charge: 0,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::C,
                atomic_number: 6,
                isotope: 13,
                charge: 0,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::H,
                atomic_number: 1,
                isotope: 0,
                charge: 1,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::O,
                atomic_number: 8,
                isotope: 0,
                charge: -2,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::Cl,
                atomic_number: 17,
                isotope: 0,
                charge: 0,
                aromatic: false,
                num_imp_h: 1,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::C,
                atomic_number: 6,
                isotope: 0,
                charge: 0,
                aromatic: false,
                num_imp_h: 1,
                chirality: Chirality::CounterClockwise,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::H,
                atomic_number: 1,
                isotope: 2,
                charge: 1,
                aromatic: false,
                num_imp_h: 0,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::C,
                atomic_number: 6,
                isotope: 13,
                charge: 0,
                aromatic: false,
                num_imp_h: 3,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::C,
                atomic_number: 6,
                isotope: 13,
                charge: 1,
                aromatic: false,
                num_imp_h: 3,
                chirality: Chirality::Undefined,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
            Atom {
                index: 0,
                atomic_symbol: AtomicSymbol::C,
                atomic_number: 6,
                isotope: 13,
                charge: 1,
                aromatic: false,
                num_imp_h: 1,
                chirality: Chirality::Clockwise,
                num_rad_electron: 0,
                ring_sizes: vec![],
                needs_update: false,
            },
        ];
        for (atom_str, atom) in atom_strs.iter().zip(atoms) {
            assert_eq!(Atom::from_str(atom_str, 0).unwrap(), atom);
        }
    }
}
