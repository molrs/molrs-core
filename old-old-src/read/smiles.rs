use std::{collections::HashMap};
use crate::mol::bond::BondType;


#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct SmilesAtom {
    pub atomic_symbol: String,
    pub aromatic: bool,
    pub isotope: u16,
    pub chirality: String,
    pub num_imp_h: u8,
    pub charge: i8,
    pub bracket: bool,
}

impl SmilesAtom {
    pub fn new(
        atomic_symbol: &str,
        aromatic: bool,
        isotope: u16,
        chirality: &str,
        num_imp_h: u8,
        charge: i8,
        bracket:bool,
    ) -> SmilesAtom {
        SmilesAtom {
            atomic_symbol: atomic_symbol.to_owned(),
            aromatic: aromatic,
            isotope: isotope,
            chirality: chirality.to_owned(),
            num_imp_h: num_imp_h,
            charge: charge,
            bracket: bracket,
        }
    }
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SmilesBond {
    pub atom_1_idx: usize,
    pub atom_2_idx: usize,
    pub bond_type: BondType,
}

impl SmilesBond {
    pub fn new(atom_1_idx: usize, atom_2_idx: usize, bond_type: BondType) -> SmilesBond {
        SmilesBond {
            atom_1_idx: atom_1_idx,
            atom_2_idx: atom_2_idx,
            bond_type: bond_type,
        }
    }
}


pub fn parse_smi(smi: &str) -> Result<(Vec<SmilesAtom>, Vec<SmilesBond>, HashMap<usize, (usize, usize)>), String> {
    let mut atoms = vec![];
    let mut bonds = vec![];
    let mut ring_closures = HashMap::new();
    if smi.chars().next().is_none() { return Ok((atoms, bonds, ring_closures)); }

    let mut atom = SmilesAtom::default();
    let mut bond_type = ' ';
    let mut bond_atom_1_idx_stack = vec![0];
    let mut c_is_in_bracket = false;
    let mut distance_from_percent_sign = 3;
    let mut double_digit_ring_closure = 0;

    let mut smi_chars = smi.chars();
    let first_char = smi_chars.next().unwrap();
    if first_char == '[' {
        c_is_in_bracket = true;
        atom.bracket = true;
    } else {
        atom.atomic_symbol.push(first_char);
        if first_char == 'b' || first_char == 'c' || first_char == 'n' ||
        first_char == 'o' || first_char == 'p' || first_char == 's' {
            atom.aromatic = true;
        } else if first_char == 'B' || first_char == 'C' || first_char == 'N' ||
        first_char == 'O' || first_char == 'P' || first_char == 'S' ||
        first_char == 'F' || first_char == 'I' || first_char == '*' {
            atom.aromatic = false;
        } else { return Err(format!("invalid first char {} in smi {}", &first_char, &smi)); }
    }

    for c in smi_chars {
        if c == '[' {
            atoms.push(atom);
            atom = SmilesAtom::default();
            atom.bracket = true;
            let len_atoms = atoms.len();
            bonds.push(SmilesBond::new(bond_atom_1_idx_stack.pop().unwrap(), len_atoms, match BondType::from_char(bond_type) {
                Some(bond_type) => bond_type,
                None => return Err(format!("invalid bond type {} in smi {}", &bond_type, &smi)),
            }));
            bond_atom_1_idx_stack.push(len_atoms);
            c_is_in_bracket = true;
        }
        if c_is_in_bracket {
            // handle bracket atoms
            if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap();
                if atom.atomic_symbol == "" { atom.isotope = atom.isotope * 10 + c_as_digit as u16; }
                else if atom.charge == 1 { atom.charge = c_as_digit as i8; }
                else if atom.charge == -1 { atom.charge = -1 * (c_as_digit as i8); }
                else if atom.num_imp_h == 1 { atom.num_imp_h = c_as_digit as u8; }
            } else if c.is_alphabetic() {
                let len_atomic_symbol = atom.atomic_symbol.chars().count();
                if len_atomic_symbol == 0 {
                    atom.atomic_symbol.push(c);
                    atom.aromatic = c.is_lowercase();
                } else if c != 'H' {
                    atom.atomic_symbol.push(c);
                } else {
                    atom.num_imp_h = 1;
                }
            } else if c == '@' {
                atom.chirality.push(c);
            } else if c == '+' {
                atom.charge = 1;
            } else if c == '-' {
                atom.charge = -1;
            } else if c == '[' || c == ']' {
            } else { return Err(format!("invalid char {} in brackets in smi {}", &c, &smi)); }
        } else {
            // handle non-bracket atoms
            if c == 'b' || c == 'c' || c == 'n' || c == 'o' || c == 'p' ||
            c == 's' {
                atoms.push(atom);
                atom = SmilesAtom::default();
                let len_atoms = atoms.len();
                bonds.push(SmilesBond::new(bond_atom_1_idx_stack.pop().unwrap(), len_atoms, match BondType::from_char(bond_type) {
                    Some(bond_type) => bond_type,
                    None => return Err(format!("invalid bond type {} in smi {}", &bond_type, &smi)),
                }));
                bond_atom_1_idx_stack.push(len_atoms);
                atom.atomic_symbol.push(c);
                atom.aromatic = true;
            } else if c == 'B' || c == 'C' || c == 'N' || c == 'O' || c == 'P' ||
            c == 'S' || c == 'F' || c == 'I' || c == '*' {
                atoms.push(atom);
                atom = SmilesAtom::default();
                let len_atoms = atoms.len();
                bonds.push(SmilesBond::new(bond_atom_1_idx_stack.pop().unwrap(), len_atoms, match BondType::from_char(bond_type) {
                    Some(bond_type) => bond_type,
                    None => return Err(format!("invalid bond type {} in smi {}", &bond_type, &smi)),
                }));
                bond_atom_1_idx_stack.push(len_atoms);
                atom.atomic_symbol.push(c);
                atom.aromatic = false;
            } else if c == 'r' && atom.atomic_symbol == "B" {
                atom.atomic_symbol.push(c);
            } else if c == 'l' && atom.atomic_symbol == "C" {
                atom.atomic_symbol.push(c);
            // handle bonds
            } else if c == '-' || c == '/' || c == '\\' || c == ':' ||
            c == '=' || c == '#' || c == '$' {
                bond_type = c;
            } else if c == '(' {
                bond_atom_1_idx_stack.push(*bond_atom_1_idx_stack.last().unwrap());
            } else if c == ')' {
                bond_atom_1_idx_stack.pop();
            // handle ring closures
            } else if c == '%' {
                distance_from_percent_sign = 0;
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap() as usize;
                if distance_from_percent_sign == 0 {
                    double_digit_ring_closure = c_as_digit * 10;
                    distance_from_percent_sign += 1;
                } else if distance_from_percent_sign == 1 {
                    double_digit_ring_closure += c_as_digit;
                    if !ring_closures.contains_key(&double_digit_ring_closure) {
                        ring_closures.insert(double_digit_ring_closure, (atoms.len(), 0));
                    } else {
                        let neighbor = ring_closures.get(&double_digit_ring_closure).unwrap().0;
                        ring_closures.insert(double_digit_ring_closure, (neighbor, atoms.len()));
                    }
                    distance_from_percent_sign += 1;
                } else {
                    if !ring_closures.contains_key(&c_as_digit) {
                        ring_closures.insert(c_as_digit, (atoms.len(), 0));
                    } else {
                        let neighbor = ring_closures.get(&c_as_digit).unwrap().0;
                        ring_closures.insert(c_as_digit, (neighbor, atoms.len()));
                    }
                }
            }
        }
        if c == ']' {
            c_is_in_bracket = false;
        }
    }
    atoms.push(atom);

    for ring_closure in ring_closures.values() {
        let bond = SmilesBond::new(ring_closure.0, ring_closure.1, BondType::Default);
        bonds.push(bond);
    }

    bonds.sort();

    Ok((atoms, bonds, ring_closures))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_smi() {
        let smi = "";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn single_non_bracket_atom() {
        let smi = "C";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "O";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("O", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "Cl";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("Cl", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn single_bracket_atom() {
        let smi = "[2H]";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("H", false, 2, "", 0, 0, true),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "[C@@]";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "@@", 0, 0, true),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "[ClH]";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("Cl", false, 0, "", 1, 0, true),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "[Cl-]";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("Cl", false, 0, "", 0, -1, true),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "[18OH-]";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("O", false, 18, "", 1, -1, true),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn linear() {
        let smi = "CC";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![
            SmilesBond::new(0, 1, BondType::Default),
        ]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "C=C";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![
            SmilesBond::new(0, 1, BondType::Double),
        ]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn parentheses() {
        let smi = "C(C)C";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![
            SmilesBond::new(0, 1, BondType::Default),
            SmilesBond::new(0, 2, BondType::Default),
        ]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "C(C)(C)C";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![
            SmilesBond::new(0, 1, BondType::Default),
            SmilesBond::new(0, 2, BondType::Default),
            SmilesBond::new(0, 3, BondType::Default),
        ]);
        assert_eq!(&ring_closures, &HashMap::new());

        let smi = "C(C(C)C)C";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![
            SmilesBond::new(0, 1, BondType::Default),
            SmilesBond::new(0, 4, BondType::Default),
            SmilesBond::new(1, 2, BondType::Default),
            SmilesBond::new(1, 3, BondType::Default),
        ]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn monocycle() {
        let smi = "C1[CH-]C1";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 1, -1, true),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![
            SmilesBond::new(0, 1, BondType::Default),
            SmilesBond::new(0, 2, BondType::Default),
            SmilesBond::new(1, 2, BondType::Default),
        ]);
        assert_eq!(&ring_closures, &HashMap::from([
            (1, (0, 2)),
        ]));

        let smi = "C%01[CH-]C1";
        let (atoms, bonds, ring_closures) = parse_smi(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
            SmilesAtom::new("C", false, 0, "", 1, -1, true),
            SmilesAtom::new("C", false, 0, "", 0, 0, false),
        ]);
        assert_eq!(&bonds, &vec![
            SmilesBond::new(0, 1, BondType::Default),
            SmilesBond::new(0, 2, BondType::Default),
            SmilesBond::new(1, 2, BondType::Default),
        ]);
        assert_eq!(&ring_closures, &HashMap::from([
            (1, (0, 2)),
        ]));
    }
}


// pub fn smi_to_atoms_and_bonds_and_ring_closures(smi: &str) -> Result<(Vec<Atom>, Vec<Bond>, HashMap<usize, (usize, usize)>), String> {
//     let mut atoms = vec![];
//     let mut bonds = vec![];
//     let mut ring_closures = HashMap::new();
    
//     // TODO: instead of an atom_str, just use a mut Atom, no more temp variables
//     let mut atom_str = String::new();
//     let mut distance_from_percent_sign = 0;
//     let mut double_digit_ring_closure = 0;
//     let mut bond_type = ' ';
//     let mut bond_start_idx_stack = vec![0];
    
//     let mut cursor_is_in_bracket = false;
//     let mut bracket_atom_symbol = String::new();
//     let mut isotope = 0;
//     let mut chirality = String::new();
//     let mut num_imp_h = 0;
//     let mut charge = 0;

//     for c in smi.chars() {
//         if c == '[' {
//             cursor_is_in_bracket = true;
//             if !atom_str.is_empty() {
//                 if atom_str.chars().next().unwrap() != '['{
//                     atoms.push(Atom::new(
//                         match AtomicSymbol::new(&atom_str) {
//                             Ok(atomic_symbol) => atomic_symbol,
//                             Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
//                         }, 0, 0, atom_str.chars().next().unwrap().is_lowercase(), 0,
//                         Chirality::Undefined, 0, atoms.len(), vec![], 0, &atom_str,
//                     ));
//                     atom_str = String::new();
//                 } else {
//                     atoms.push(Atom::new(
//                         match AtomicSymbol::new(&bracket_atom_symbol) {
//                             Ok(atomic_symbol) => atomic_symbol,
//                             Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
//                         }, isotope, charge, bracket_atom_symbol.chars().next().unwrap().is_lowercase(), num_imp_h,
//                         match Chirality::new(&chirality) {
//                             Ok(chirality) => chirality,
//                             Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
//                         }, 0, atoms.len(), vec![], 0, &atom_str,
//                     ));
//                     atom_str = String::new();
//                     bracket_atom_symbol = String::new();
//                     isotope = 0;
//                     charge = 0;
//                     num_imp_h = 0;
//                     chirality = String::new();
//                 }
//                 let num_atoms = atoms.len();
//                 bonds.push(Bond::new(bond_start_idx_stack.pop().unwrap(), num_atoms, match BondType::new(bond_type) {
//                     Ok(bond_type) => bond_type,
//                     Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
//                 }));
//                 bond_type = ' ';
//                 bond_start_idx_stack.push(num_atoms);
//             }
//         }
//         if cursor_is_in_bracket {
//             atom_str.push(c);
//             if c.is_numeric() {
//                 let c_as_digit = c.to_digit(10).unwrap() as usize;
//                 if bracket_atom_symbol.is_empty() {
//                     isotope = (isotope * 10) + c_as_digit;
//                 } else if charge == -1 {
//                     charge = c_as_digit as isize * -1;
//                 } else if charge == 1 {
//                     charge = c_as_digit as isize;
//                 } else if num_imp_h == 1 {
//                     num_imp_h = c_as_digit;
//                 } else { return Err(format!("unexpected numeric in bracket atom in smi {}", &smi)) }
//             } else if c.is_alphabetic() {
//                 let len_bracket_atom_symbol = bracket_atom_symbol.chars().count();
//                 if len_bracket_atom_symbol == 0 {
//                     bracket_atom_symbol.push(c);
//                 } else if len_bracket_atom_symbol == 1 && c != 'H' {
//                     bracket_atom_symbol.push(c);
//                 } else if c == 'H' {
//                     num_imp_h = 1;
//                 }
//             } else if c == '@' {
//                 chirality.push(c);
//             } else if c == '-' {
//                 charge = -1;
//             } else if c == '+' {
//                 charge = 1;
//             } else if c == '[' || c == ']' {
//             } else { return Err(format!("invalid char {} in bracket atom in smi {}", &c, &smi)) }
//         } else {
//             if c == '(' {
//                 bond_start_idx_stack.push(*bond_start_idx_stack.last().unwrap());
//             } else if c == ')' { bond_start_idx_stack.pop();
//             } else if c == 'b' || c == 'B' || c == 'c' || c == 'C' || c == 'n' || c == 'N' ||
//             c == 'o' || c == 'O' || c == 'p' || c == 'P' || c == 's' || c == 'S' ||
//             c == 'F' || c == 'I' || c == '*' {
//                 if !atom_str.is_empty() {
//                     if atom_str.chars().next().unwrap() != '['{
//                         atoms.push(Atom::new(
//                             match AtomicSymbol::new(&atom_str) {
//                                 Ok(atomic_symbol) => atomic_symbol,
//                                 Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
//                             }, 0, 0, atom_str.chars().next().unwrap().is_lowercase(), 0,
//                             Chirality::Undefined, 0, atoms.len(), vec![], 0, &atom_str,
//                         ));
//                         atom_str = String::new();
//                     } else {
//                         atoms.push(Atom::new(
//                             match AtomicSymbol::new(&bracket_atom_symbol) {
//                                 Ok(atomic_symbol) => atomic_symbol,
//                                 Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
//                             }, isotope, charge, bracket_atom_symbol.chars().next().unwrap().is_lowercase(), num_imp_h,
//                             match Chirality::new(&chirality) {
//                                 Ok(chirality) => chirality,
//                                 Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
//                             }, 0, atoms.len(), vec![], 0, &atom_str,
//                         ));
//                         atom_str = String::new();
//                         bracket_atom_symbol = String::new();
//                         isotope = 0;
//                         charge = 0;
//                         num_imp_h = 0;
//                         chirality = String::new();
//                     }
//                     let num_atoms = atoms.len();
//                     bonds.push(Bond::new(bond_start_idx_stack.pop().unwrap(), num_atoms, match BondType::new(bond_type) {
//                         Ok(bond_type) => bond_type,
//                         Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
//                     }));
//                     bond_type = ' ';
//                     bond_start_idx_stack.push(num_atoms);
//                 }
//                 atom_str.push(c);
//             } else if c == 'r' || c == 'l' {
//                 atom_str.push(c);
//             } else if c == '%' {
//                 distance_from_percent_sign = 2;
//             } else if c.is_numeric() {
//                 let c_as_digit = c.to_digit(10).unwrap() as usize;
//                 if distance_from_percent_sign == 2 {
//                     distance_from_percent_sign -= 1;
//                     double_digit_ring_closure += c_as_digit * 10;
//                 } else if distance_from_percent_sign == 1 {
//                     distance_from_percent_sign -= 1;
//                     double_digit_ring_closure += c_as_digit;
//                     if !ring_closures.contains_key(&double_digit_ring_closure) {
//                         ring_closures.insert(double_digit_ring_closure, (atoms.len(), 0));
//                     } else {
//                         let neighbor = ring_closures.get(&double_digit_ring_closure).unwrap().0;
//                         ring_closures.insert(double_digit_ring_closure, (neighbor, atoms.len()));
//                     }
//                     double_digit_ring_closure = 0;
//                 } else if distance_from_percent_sign == 0 {
//                     if !ring_closures.contains_key(&c_as_digit) {
//                         ring_closures.insert(c_as_digit, (atoms.len(), 0));
//                     } else {
//                         let neighbor = ring_closures.get(&c_as_digit).unwrap().0;
//                         ring_closures.insert(c_as_digit, (neighbor, atoms.len()));
//                     }
//                 }
//             } else if c == '-' || c == '/' || c == '\\' || c == ':' || c == '=' ||
//             c == '#' || c == '$' || c == '~' {
//                 bond_type = c;
//             } else { return Err(format!("invalid char {} in smi {}", &c, &smi)); }
//         }
//         if c == ']' {
//             cursor_is_in_bracket = false;
//         }
//     }
//     if !atom_str.is_empty() {
//         if atom_str.chars().next().unwrap() != '['{
//             atoms.push(Atom::new(
//                 match AtomicSymbol::new(&atom_str) {
//                     Ok(atomic_symbol) => atomic_symbol,
//                     Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
//                 }, 0, 0, atom_str.chars().next().unwrap().is_lowercase(), 0,
//                 Chirality::Undefined, 0, atoms.len(), vec![], 0, &atom_str,
//             ));
//         } else {
//             atoms.push(Atom::new(
//                 match AtomicSymbol::new(&bracket_atom_symbol) {
//                     Ok(atomic_symbol) => atomic_symbol,
//                     Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
//                 }, isotope, charge, bracket_atom_symbol.chars().next().unwrap().is_lowercase(), num_imp_h,
//                 match Chirality::new(&chirality) {
//                     Ok(chirality) => chirality,
//                     Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
//                 }, 0, atoms.len(), vec![], 0, &atom_str,
//             ));
//         }
//     }

//     for ring_closure in ring_closures.values() {
//         let i = ring_closure.0;
//         let j = ring_closure.1;
//         let bond = Bond::from_char(i, j, ' ').unwrap();
//         bonds.push(bond);
//     }
//     bonds.sort();

//     for (i, atom) in atoms.iter_mut().enumerate() {
//         atom.neighbor_idxs = get_neighbor_idxs(i, &bonds);
//     }

//     Ok((atoms, bonds, ring_closures))
// }

// /// Gives the indices of the neighbors of an atom specified by its atom_idx.
// /// 
// /// # Arguments
// /// * `atom_idx`: The index of the atom in the molecule.
// /// * `bonds`: A vector of all the bonds in a molecule.
// /// 
// /// # Examples
// /// ```
// /// let smi = "CC";
// /// let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds(&smi).unwrap();
// /// let neighbor_idxs = get_neighbor_idxs(0, &bonds);
// /// assert_eq!(&neighbor_idxs, &vec![1]);
// /// ```
// fn get_neighbor_idxs(atom_idx: usize, bonds: &Vec<Bond>) -> Vec<usize> {
//     let mut neighbor_idxs = vec![];

//     for bond in bonds.iter() {
//         if bond.atom_idx_1 == atom_idx { neighbor_idxs.push(bond.atom_idx_2); }
//         if bond.atom_idx_2 == atom_idx { neighbor_idxs.push(bond.atom_idx_1); }
//     }

//     neighbor_idxs
// }


// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn parse_empty_smi() {
//         let smi = "";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![]);
//         assert_eq!(&bonds, &vec![]);
//         assert_eq!(&ring_closures, &HashMap::new());
//     }

//     #[test]
//     fn parse_linear_smi() {
//         let smi = "C";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![]);
//         assert_eq!(&ring_closures, &HashMap::new());
//         let smi = "ClCCl";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::Cl, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1], 0, "Cl"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0, 2], 0, "C"),
//             Atom::new(AtomicSymbol::Cl, 0, 0, false, 0, Chirality::Undefined, 0, 2, vec![1], 0, "Cl"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Default),
//             Bond::new(1, 2, BondType::Default),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//         let smi = "OCCOCCOCCO";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::O, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1], 0, "O"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0, 2], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 2, vec![1, 3], 0, "C"),
//             Atom::new(AtomicSymbol::O, 0, 0, false, 0, Chirality::Undefined, 0, 3, vec![2, 4], 0, "O"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 4, vec![3, 5], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 5, vec![4, 6], 0, "C"),
//             Atom::new(AtomicSymbol::O, 0, 0, false, 0, Chirality::Undefined, 0, 6, vec![5, 7], 0, "O"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 7, vec![6, 8], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 8, vec![7, 9], 0, "C"),
//             Atom::new(AtomicSymbol::O, 0, 0, false, 0, Chirality::Undefined, 0, 9, vec![8], 0, "O"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Default),
//             Bond::new(1, 2, BondType::Default),
//             Bond::new(2, 3, BondType::Default),
//             Bond::new(3, 4, BondType::Default),
//             Bond::new(4, 5, BondType::Default),
//             Bond::new(5, 6, BondType::Default),
//             Bond::new(6, 7, BondType::Default),
//             Bond::new(7, 8, BondType::Default),
//             Bond::new(8, 9, BondType::Default),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//     }

//     #[test]
//     fn parse_bond_smi() {
//         let smi = "C-C";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Single),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//         let smi = "C=C";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Double),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//         let smi = "C#C";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Triple),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//     }

//     #[test]
//     fn parse_bracket_atom() {
//         let smi = "[13CH3+]";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 13, 1, false, 3, Chirality::Undefined, 0, 0, vec![], 0, "[13CH3+]"),
//         ]);
//         assert_eq!(&bonds, &vec![]);
//         assert_eq!(&ring_closures, &HashMap::new());
//     }

//     #[test]
//     fn parse_parentheses() {
//         let smi = "C(C)C";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1, 2], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 2, vec![0], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Default),
//             Bond::new(0, 2, BondType::Default),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//         let smi = "C(C(C)C)C";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1, 4], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0, 2, 3], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 2, vec![1], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 3, vec![1], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 4, vec![0], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Default),
//             Bond::new(0, 4, BondType::Default),
//             Bond::new(1, 2, BondType::Default),
//             Bond::new(1, 3, BondType::Default),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//         let smi = "C(C)(C)C";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1, 2, 3], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 2, vec![0], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 3, vec![0], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Default),
//             Bond::new(0, 2, BondType::Default),
//             Bond::new(0, 3, BondType::Default),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::new());
//     }

//     #[test]
//     fn parse_ring_bonds() {
//         let smi = "C1C(C2)C12";
//         let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
//         assert_eq!(&atoms, &vec![
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1, 3], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![0, 2, 3], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 2, vec![1, 3], 0, "C"),
//             Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 3, vec![0, 1, 2], 0, "C"),
//         ]);
//         assert_eq!(&bonds, &vec![
//             Bond::new(0, 1, BondType::Default),
//             Bond::new(0, 3, BondType::Default),
//             Bond::new(1, 2, BondType::Default),
//             Bond::new(1, 3, BondType::Default),
//             Bond::new(2, 3, BondType::Default),
//         ]);
//         assert_eq!(&ring_closures, &HashMap::from([
//             (1, (0, 3)),
//             (2, (2, 3)),
//         ]));
//     }

//     #[test]
//     fn test_get_neighbor_idxs() {
//         let bonds = vec![
//             Bond::new(0, 1, BondType::Single),
//             Bond::new(0, 2, BondType::Single),
//             Bond::new(0, 3, BondType::Double),
//         ];
//         let neighbor_idxs = get_neighbor_idxs(0, &bonds);

//         assert_eq!(neighbor_idxs, vec![1, 2, 3]);
//     }
// }