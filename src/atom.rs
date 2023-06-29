use std::str::FromStr;

use crate::element::{Element, ElementParseError};
pub struct AtomParseError {
    pub details: String,
}

impl From<ElementParseError> for AtomParseError {
    fn from(element_parse_error: ElementParseError) -> Self {
        AtomParseError {
            details: element_parse_error.details,
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum PointChirality {
    #[default]
    Undefined,
    Clockwise,
    CounterClockwise,
}

impl FromStr for PointChirality {
    type Err = AtomParseError;

    fn from_str(point_chirality: &str) -> Result<Self, Self::Err> {
        match point_chirality {
            "" => Ok(PointChirality::Undefined),
            "@" => Ok(PointChirality::CounterClockwise),
            "@@" => Ok(PointChirality::Clockwise),
            _ => Err(AtomParseError {
                details: format!("invalid point chirality {}", point_chirality),
            }),
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Atom {
    pub element: Element,
    pub isotope: u16,
    pub charge: i8,
    pub aromatic: bool,
    pub num_implicit_hydrogens: u8,
    pub num_radical_electrons: u8,
    pub point_chirality: PointChirality,
}

impl FromStr for Atom {
    type Err = AtomParseError;

    fn from_str(atom_str: &str) -> Result<Self, Self::Err> {
        if atom_str.starts_with('[') {
            atom_from_str_bracket(atom_str)
        } else {
            atom_from_str_non_bracket(atom_str)
        }
    }
}

enum NextValueInBracket {
    Isotope,
    Charge(bool),
    NumImpH,
}

fn atom_from_str_bracket(atom_str: &str) -> Result<Atom, AtomParseError> {
    let mut atomic_symbol = String::new();
    let mut atom = Atom::default();
    let mut next_value = NextValueInBracket::Isotope;
    for c in atom_str.chars() {
        if c == '[' {
        } else if c == ']' {
            break;
        } else if c.is_numeric() {
            let c_as_int = c.to_digit(10).unwrap();
            match next_value {
                NextValueInBracket::Isotope => {
                    atom.isotope = atom.isotope * 10 + c_as_int as u16;
                }
                NextValueInBracket::Charge(is_positive) => {
                    if is_positive {
                        atom.charge = c_as_int as i8;
                    } else {
                        atom.charge = -(c_as_int as i8);
                    }
                }
                NextValueInBracket::NumImpH => {
                    atom.num_implicit_hydrogens = c_as_int as u8;
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
                atom.num_implicit_hydrogens = 1;
                next_value = NextValueInBracket::NumImpH;
            } else {
            }
        } else if c == '@' {
            atom.point_chirality = match atom.point_chirality {
                PointChirality::Undefined => PointChirality::CounterClockwise,
                PointChirality::CounterClockwise => PointChirality::Clockwise,
                _ => {
                    return Err(AtomParseError {
                        details: format!("problem with point chirality in {}", atom_str),
                    })
                }
            }
        } else if c == '+' {
            atom.charge = 1;
            next_value = NextValueInBracket::Charge(true);
        } else if c == '-' {
            atom.charge = -1;
            next_value = NextValueInBracket::Charge(false);
        } else {
        }
    }
    atom.element = Element::from_str(&atomic_symbol)?;

    Ok(atom)
}

fn atom_from_str_non_bracket(atom_str: &str) -> Result<Atom, AtomParseError> {
    let mut atomic_symbol = String::new();
    let mut atom = Atom::default();
    for c in atom_str.chars() {
        if !c.is_alphabetic() && c != '*' {
            return Err(AtomParseError {
                details: format!("invalid non-bracket char {}", c),
            });
        }
        if atomic_symbol.is_empty() {
            atom.aromatic = c.is_lowercase();
        }
        atomic_symbol.push(c);
    }
    atom.element = Element::from_str(&atomic_symbol)?;

    Ok(atom)
}
