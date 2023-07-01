use std::{collections::HashMap, str::FromStr};

use crate::{atom::Atom, bond::Bond};

#[derive(Debug)]
pub struct SmilesParseError {
    pub details: String,
}

#[derive(Debug)]
enum RingClosure {
    OneAtom(usize),
    Used,
}

/// Struct to handle parsing of SMILES (and related formats) to atoms and bonds.
///
/// ```
/// use molrs::{smiles::SmilesParser, atom::Atom, bond::Bond};
///
/// let smiles_parser: SmilesParser = "C-C".parse().unwrap();
/// let atoms: Vec<Atom> = smiles_parser.atoms().unwrap();
/// let bonds: Vec<Bond> = smiles_parser.bonds().unwrap();
/// ```
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SmilesParser {
    pub smi: String,
    pub atom_strs: Vec<String>,
    pub bond_strs: Vec<String>,
    pub ring_closure_strs: Vec<String>,
}

impl FromStr for SmilesParser {
    type Err = SmilesParseError;

    /// Parse a SMILES (and related formats) to atom_strs, bond_strs, and
    /// ring_closure_strs.
    ///
    /// An atom_str represents a single atom. "c" or "Cl" can be an atom_str.
    /// [14C-H1] can also be an atom str. A bond_str represents a bond. An empty
    /// bond_str indicates a default bond, "-" means single bond, "=" means
    /// double bond, and so on. A ring_closure_str represents additional bonds
    /// needed to form ring closures. These are typically represented as single-
    /// digit numbers or %-prepended double-digit numbers that follow an atom.
    ///
    /// ```
    /// use molrs::{smiles::SmilesParser, atom::Atom, bond::Bond};
    ///
    /// let smiles_parser: SmilesParser = "C-C".parse().unwrap();
    ///
    /// let atom_strs = ["C", "C"];
    /// let atom_strs = atom_strs.into_iter().map(|s| s.to_owned()).collect::<Vec<String>>();
    /// assert_eq!(smiles_parser.atom_strs, atom_strs);
    ///
    /// let bond_strs = ["-"];
    /// let bond_strs = bond_strs.into_iter().map(|s| s.to_owned()).collect::<Vec<String>>();
    /// assert_eq!(smiles_parser.bond_strs, bond_strs)
    /// ```
    fn from_str(smi: &str) -> Result<Self, Self::Err> {
        let mut atom_strs: Vec<String> = vec![];
        let mut bond_strs = vec![];
        let mut ring_closure_strs = vec![];

        let mut c_is_in_bracket = false;
        for c in smi.chars() {
            if c_is_in_bracket {
                atom_strs.last_mut().unwrap().push(c);
                if c == ']' {
                    c_is_in_bracket = false;
                }
            } else if c == '[' {
                c_is_in_bracket = true;
                atom_strs.push(String::from(c));
                bond_strs.push(String::new());
            } else if c == '-'
                || c == '/'
                || c == '\\'
                || c == ':'
                || c == '='
                || c == '#'
                || c == '$'
                || c == '('
                || c == ')'
                || c == '%'
                || c == '.'
                || c.is_numeric()
            {
                bond_strs.last_mut().unwrap().push(c);
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
                atom_strs.push(String::from(c));
                bond_strs.push(String::new());
            } else if c == 'l' || c == 'r' {
                atom_strs.last_mut().unwrap().push(c);
            } else {
                return Err(SmilesParseError {
                    details: format!("{} | invalid char {}", smi, c),
                });
            }
        }

        for bond_str in bond_strs.iter_mut() {
            let mut split_index = bond_str.chars().count();
            for c in bond_str.chars().rev() {
                if c.is_numeric() {
                    break;
                }
                split_index -= 1;
            }
            ring_closure_strs.push(bond_str[..split_index].to_owned());
            *bond_str = bond_str[split_index..].to_owned();
        }
        bond_strs.pop();

        Ok(SmilesParser {
            smi: smi.to_owned(),
            atom_strs,
            bond_strs,
            ring_closure_strs,
        })
    }
}

impl SmilesParser {
    /// Parses atom_strs: Vec<String> to atoms: Vec<Atom>.
    ///
    /// ```
    /// use molrs::{smiles::SmilesParser, atom::Atom, element::Element};
    ///
    /// let smiles_parser: SmilesParser = "C".parse().unwrap();
    ///
    /// let mut atom = Atom::default();
    /// atom.element = Element::C;
    /// let atoms = vec![atom];
    /// assert_eq!(smiles_parser.atoms().unwrap(), atoms)
    /// ```
    pub fn atoms(&self) -> Result<Vec<Atom>, SmilesParseError> {
        match self.atom_strs.iter().map(|s| s.parse()).try_collect() {
            Ok(atoms) => Ok(atoms),
            Err(atom_parse_error) => Err(SmilesParseError {
                details: format!("{} | {}", self.smi, atom_parse_error.details),
            }),
        }
    }

    /// Parses bond_strs: Vec<String> and ring_closure_strs: Vec<String> to
    /// bonds: Vec<Bond>.
    ///
    /// Handles both bond_strs for linear/branched bonds and ring_closure_strs
    /// for ring-forming bonds. The first phase uses a stack to keep track of
    /// which atom index is the current root of the tree. The second phase uses
    /// a hash map to map a ring index to the atom index that needs to be
    /// connected.
    ///
    /// ```
    /// use molrs::{smiles::SmilesParser, bond::Bond};
    ///
    /// let smiles_parser: SmilesParser = "C1CC1".parse().unwrap();
    ///
    /// let bonds = vec![
    ///     Bond::new(0, 1, ' '),
    ///     Bond::new(1, 2, ' '),
    ///     Bond::new(0, 2, ' '),
    /// ];
    /// let bonds: Vec<Bond> = bonds.into_iter().map(|res| res.unwrap()).collect();
    /// assert_eq!(smiles_parser.bonds().unwrap(), bonds);
    /// ```
    pub fn bonds(&self) -> Result<Vec<Bond>, SmilesParseError> {
        let mut bonds = vec![];

        let mut source_atom_index_stack = vec![0];
        for (i, bond_str) in self.bond_strs.iter().enumerate() {
            let target_atom_index = i + 1;
            let mut bond_char = ' ';
            for c in bond_str.chars() {
                if c == '(' {
                    source_atom_index_stack.push(*source_atom_index_stack.last().unwrap());
                } else if c == ')' {
                    source_atom_index_stack.pop();
                } else if c == '-'
                    || c == '/'
                    || c == '\\'
                    || c == ':'
                    || c == '='
                    || c == '#'
                    || c == '$'
                    || c == '.'
                {
                    bond_char = c;
                } else {
                    return Err(SmilesParseError {
                        details: format!("{} | invalid char {}", self.smi, c),
                    });
                }
            }
            if bond_char != '.' {
                bonds.push(
                    match Bond::new(
                        source_atom_index_stack.pop().unwrap(),
                        target_atom_index,
                        bond_char,
                    ) {
                        Ok(bond) => bond,
                        Err(bond_parse_error) => {
                            return Err(SmilesParseError {
                                details: format!("{} | {}", self.smi, bond_parse_error.details),
                            })
                        }
                    },
                );
                source_atom_index_stack.push(target_atom_index);
            }
        }

        let mut ring_closures = HashMap::new();
        for (i, ring_closure_str) in self.ring_closure_strs.iter().enumerate() {
            let source_atom_index = i;
            let mut ring_indices = vec![];
            let mut bond_char = ' ';
            let mut distance_from_percent = 2;
            for c in ring_closure_str.chars() {
                if distance_from_percent == 0 {
                    distance_from_percent += 1;
                    ring_indices.push(c.to_digit(10).unwrap());
                } else if distance_from_percent == 1 {
                    distance_from_percent += 1;
                    let ring_index = ring_indices.last_mut().unwrap();
                    *ring_index *= 10;
                    *ring_index += c.to_digit(10).unwrap();
                } else if c == '-'
                    || c == '/'
                    || c == '\\'
                    || c == ':'
                    || c == '='
                    || c == '#'
                    || c == '$'
                {
                    bond_char = c;
                } else if c == '%' {
                    distance_from_percent = 0;
                } else if c.is_numeric() {
                    ring_indices.push(c.to_digit(10).unwrap());
                } else {
                    return Err(SmilesParseError {
                        details: format!("{} | invalid char {}", self.smi, c),
                    });
                }
            }
            if ring_indices.is_empty() {
                continue;
            }
            for ring_index in ring_indices {
                if let Some(ring_closure) = ring_closures.get(&ring_index) {
                    match ring_closure {
                        RingClosure::OneAtom(target_atom_index) => {
                            bonds.push(
                                match Bond::new(*target_atom_index, source_atom_index, bond_char) {
                                    Ok(bond) => bond,
                                    Err(bond_parse_error) => {
                                        return Err(SmilesParseError {
                                            details: format!(
                                                "{} | {}",
                                                self.smi, bond_parse_error.details
                                            ),
                                        })
                                    }
                                },
                            );
                            ring_closures.insert(ring_index, RingClosure::Used);
                        }
                        RingClosure::Used => {
                            return Err(SmilesParseError {
                                details: format!(
                                    "{} | repeated ring index {}",
                                    self.smi, ring_index
                                ),
                            });
                        }
                    }
                } else {
                    ring_closures.insert(ring_index, RingClosure::OneAtom(source_atom_index));
                }
            }
        }

        for ring_index in ring_closures.keys() {
            let ring_closure = ring_closures.get(ring_index).unwrap();
            match ring_closure {
                RingClosure::OneAtom(_) => {
                    return Err(SmilesParseError {
                        details: format!("{} | unused ring index {}", self.smi, ring_index),
                    })
                }
                RingClosure::Used => (),
            }
        }

        Ok(bonds)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smiles_parser_from_str() {
        let data = [
            ("C", vec!["C"], vec![], vec![""]),
            ("CC", vec!["C", "C"], vec![""], vec!["", ""]),
            ("C=C", vec!["C", "C"], vec!["="], vec!["", ""]),
            ("CC=O", vec!["C", "C", "O"], vec!["", "="], vec!["", "", ""]),
            ("C[H]", vec!["C", "[H]"], vec![""], vec!["", ""]),
            ("[18O]C", vec!["[18O]", "C"], vec![""], vec!["", ""]),
            (
                "C1CC=1",
                vec!["C", "C", "C"],
                vec!["", ""],
                vec!["1", "", "=1"],
            ),
            (
                "C%10C[N-]%10",
                vec!["C", "C", "[N-]"],
                vec!["", ""],
                vec!["%10", "", "%10"],
            ),
            (
                "C(F)C",
                vec!["C", "F", "C"],
                vec!["(", ")"],
                vec!["", "", ""],
            ),
            (
                "C(C)(C)C",
                vec!["C", "C", "C", "C"],
                vec!["(", ")(", ")"],
                vec!["", "", "", ""],
            ),
            (
                "C(C(C)C)C",
                vec!["C", "C", "C", "C", "C"],
                vec!["(", "(", ")", ")"],
                vec!["", "", "", "", ""],
            ),
        ];
        for (smi, atom_strs, bond_strs, ring_closure_strs) in data {
            let atom_strs: Vec<String> = atom_strs.into_iter().map(|s| s.to_owned()).collect();
            let bond_strs: Vec<String> = bond_strs.into_iter().map(|s| s.to_owned()).collect();
            let ring_closure_strs: Vec<String> = ring_closure_strs
                .into_iter()
                .map(|s| s.to_owned())
                .collect();
            assert_eq!(
                SmilesParser::from_str(smi).unwrap(),
                SmilesParser {
                    smi: smi.to_owned(),
                    atom_strs,
                    bond_strs,
                    ring_closure_strs
                }
            );
        }
    }
}
