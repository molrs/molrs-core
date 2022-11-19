use std::{collections::HashMap};
use crate::{periodic_table::{AtomicSymbol}, mol::{atom::{Atom, Chirality}, bond::{Bond, BondType}}};


pub fn smi_to_atoms_and_bonds_and_ring_closures(smi: &str) -> Result<(Vec<Atom>, Vec<Bond>, HashMap<usize, (usize, usize)>), String> {
    let mut atoms = vec![];
    let mut bonds = vec![];
    let mut ring_closures = HashMap::new();
    
    let mut atom_str = String::new();
    let mut distance_from_percent_sign = 0;
    let mut double_digit_ring_closure = 0;
    let mut bond_type = ' ';
    let mut bond_start_idx_stack = vec![0];
    
    let mut cursor_is_in_bracket = false;
    let mut bracket_atom_symbol = String::new();
    let mut isotope = 0;
    let mut chirality = String::new();
    let mut num_imp_h = 0;
    let mut charge = 0;

    for c in smi.chars() {
        if c == '[' {
            cursor_is_in_bracket = true;
            if !atom_str.is_empty() {
                if atom_str.chars().next().unwrap() != '['{
                    atoms.push(Atom::new(
                        match AtomicSymbol::new(&atom_str) {
                            Ok(atomic_symbol) => atomic_symbol,
                            Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
                        }, 0, 0, atom_str.chars().next().unwrap().is_lowercase(), 0,
                        Chirality::Undefined, 0, atoms.len(), vec![], 0, &atom_str,
                    ));
                    atom_str = String::new();
                } else {
                    atoms.push(Atom::new(
                        match AtomicSymbol::new(&bracket_atom_symbol) {
                            Ok(atomic_symbol) => atomic_symbol,
                            Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
                        }, isotope, charge, bracket_atom_symbol.chars().next().unwrap().is_lowercase(), num_imp_h,
                        match Chirality::new(&chirality) {
                            Ok(chirality) => chirality,
                            Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
                        }, 0, atoms.len(), vec![], 0, &atom_str,
                    ));
                    atom_str = String::new();
                    bracket_atom_symbol = String::new();
                    isotope = 0;
                    charge = 0;
                    num_imp_h = 0;
                    chirality = String::new();
                }
                let num_atoms = atoms.len();
                bonds.push(Bond::new(bond_start_idx_stack.pop().unwrap(), num_atoms, match BondType::new(bond_type) {
                    Ok(bond_type) => bond_type,
                    Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
                }));
                bond_type = ' ';
                bond_start_idx_stack.push(num_atoms);
            }
        }
        if cursor_is_in_bracket {
            atom_str.push(c);
            if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap() as usize;
                if bracket_atom_symbol.is_empty() {
                    isotope = (isotope * 10) + c_as_digit;
                } else if charge == -1 {
                    charge = c_as_digit as isize * -1;
                } else if charge == 1 {
                    charge = c_as_digit as isize;
                } else if num_imp_h == 1 {
                    num_imp_h = c_as_digit;
                } else { return Err(format!("unexpected numeric in bracket atom in smi {}", &smi)) }
            } else if c.is_alphabetic() {
                let len_bracket_atom_symbol = bracket_atom_symbol.chars().count();
                if len_bracket_atom_symbol == 0 {
                    bracket_atom_symbol.push(c);
                } else if len_bracket_atom_symbol == 1 && c != 'H' {
                    bracket_atom_symbol.push(c);
                } else if c == 'H' {
                    num_imp_h = 1;
                }
            } else if c == '@' {
                chirality.push(c);
            } else if c == '-' {
                charge = -1;
            } else if c == '+' {
                charge = 1;
            } else if c == '[' || c == ']' {
            } else { return Err(format!("invalid char {} in bracket atom in smi {}", &c, &smi)) }
        } else {
            if c == '(' {
                bond_start_idx_stack.push(*bond_start_idx_stack.last().unwrap());
            } else if c == ')' { bond_start_idx_stack.pop();
            } else if c == 'b' || c == 'B' || c == 'c' || c == 'C' || c == 'n' || c == 'N' ||
            c == 'o' || c == 'O' || c == 'p' || c == 'P' || c == 's' || c == 'S' ||
            c == 'F' || c == 'I' || c == '*' {
                if !atom_str.is_empty() {
                    if atom_str.chars().next().unwrap() != '['{
                        atoms.push(Atom::new(
                            match AtomicSymbol::new(&atom_str) {
                                Ok(atomic_symbol) => atomic_symbol,
                                Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
                            }, 0, 0, atom_str.chars().next().unwrap().is_lowercase(), 0,
                            Chirality::Undefined, 0, atoms.len(), vec![], 0, &atom_str,
                        ));
                        atom_str = String::new();
                    } else {
                        atoms.push(Atom::new(
                            match AtomicSymbol::new(&bracket_atom_symbol) {
                                Ok(atomic_symbol) => atomic_symbol,
                                Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
                            }, isotope, charge, bracket_atom_symbol.chars().next().unwrap().is_lowercase(), num_imp_h,
                            match Chirality::new(&chirality) {
                                Ok(chirality) => chirality,
                                Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
                            }, 0, atoms.len(), vec![], 0, &atom_str,
                        ));
                        atom_str = String::new();
                        bracket_atom_symbol = String::new();
                        isotope = 0;
                        charge = 0;
                        num_imp_h = 0;
                        chirality = String::new();
                    }
                    let num_atoms = atoms.len();
                    bonds.push(Bond::new(bond_start_idx_stack.pop().unwrap(), num_atoms, match BondType::new(bond_type) {
                        Ok(bond_type) => bond_type,
                        Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
                    }));
                    bond_type = ' ';
                    bond_start_idx_stack.push(num_atoms);
                }
                atom_str.push(c);
            } else if c == 'r' || c == 'l' {
                atom_str.push(c);
            } else if c == '%' {
                distance_from_percent_sign = 2;
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap() as usize;
                if distance_from_percent_sign == 2 {
                    distance_from_percent_sign -= 1;
                    double_digit_ring_closure += c_as_digit * 10;
                } else if distance_from_percent_sign == 1 {
                    distance_from_percent_sign -= 1;
                    double_digit_ring_closure += c_as_digit;
                    if !ring_closures.contains_key(&double_digit_ring_closure) {
                        ring_closures.insert(double_digit_ring_closure, (atoms.len(), 0));
                    } else {
                        let neighbor = ring_closures.get(&double_digit_ring_closure).unwrap().0;
                        ring_closures.insert(double_digit_ring_closure, (neighbor, atoms.len()));
                    }
                    double_digit_ring_closure = 0;
                } else if distance_from_percent_sign == 0 {
                    if !ring_closures.contains_key(&c_as_digit) {
                        ring_closures.insert(c_as_digit, (atoms.len(), 0));
                    } else {
                        let neighbor = ring_closures.get(&c_as_digit).unwrap().0;
                        ring_closures.insert(c_as_digit, (neighbor, atoms.len()));
                    }
                }
            } else if c == '-' || c == '/' || c == '\\' || c == ':' || c == '=' ||
            c == '#' || c == '$' || c == '~' {
                bond_type = c;
            } else { return Err(format!("invalid char {} in smi {}", &c, &smi)); }
        }
        if c == ']' {
            cursor_is_in_bracket = false;
        }
    }
    if !atom_str.is_empty() {
        if atom_str.chars().next().unwrap() != '['{
            atoms.push(Atom::new(
                match AtomicSymbol::new(&atom_str) {
                    Ok(atomic_symbol) => atomic_symbol,
                    Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
                }, 0, 0, atom_str.chars().next().unwrap().is_lowercase(), 0,
                Chirality::Undefined, 0, atoms.len(), vec![], 0, &atom_str,
            ));
        } else {
            atoms.push(Atom::new(
                match AtomicSymbol::new(&bracket_atom_symbol) {
                    Ok(atomic_symbol) => atomic_symbol,
                    Err(problem) => return Err(format!("problem: {}, in smi {}", &problem, &smi)),
                }, isotope, charge, bracket_atom_symbol.chars().next().unwrap().is_lowercase(), num_imp_h,
                match Chirality::new(&chirality) {
                    Ok(chirality) => chirality,
                    Err(problem) => return Err(format!("problem {} in smi {}", &problem, &smi)),
                }, 0, atoms.len(), vec![], 0, &atom_str,
            ));
        }
    }

    for ring_closure in ring_closures.values() {
        let i = ring_closure.0;
        let j = ring_closure.1;
        let bond = Bond::from_char(i, j, ' ').unwrap();
        bonds.push(bond);
    }
    bonds.sort();

    // everything below this line should be moved to Mol normalization

    perceive_default_bonds(&atoms, &mut bonds);

    for (i, atom) in atoms.iter_mut().enumerate() {
        let (valence, neighbor_idxs) = get_explicit_valence_and_neighbor_idxs(i, &bonds);
        // if !atom.aromatic && atom.as_smi.chars().next().unwrap() != '[' {
        if !atom.aromatic {
            let max_allowed_valence = atom.max_allowed_valence();
            if max_allowed_valence.is_some() {
                let max_allowed_valence = max_allowed_valence.unwrap();
                if valence > max_allowed_valence { return Err(format!("atom exceeded maximum valence in smi {}, {}", &smi, &i)); }
                let default_num_imp_h = atom.atomic_symbol.num_imp_h();
                if valence > default_num_imp_h {
                    atom.num_imp_h = (max_allowed_valence - valence) % 2;
                } else {
                    atom.num_imp_h = default_num_imp_h - valence;
                }
            } else {
                return Err(format!("unimplemented atomic state {} in smi {}", &atom.as_smi, &smi))
                // need to think harder about how to handle organometallic complexes
                // unimplemented!();
            }
        }
        atom.neighbor_idxs = neighbor_idxs;
    }

    Ok((atoms, bonds, ring_closures))
}


/// Converts default BondTypes to single or aromatic in-place.
/// 
/// # Arguments
/// * `atoms`: A vector of all the atoms in a molecule.
/// * `bonds`: A vector of all the bonds in a molecule.
/// 
/// # Examples
/// ```
/// let smi = "Cc";
/// let (atoms, mut bonds, ring_closures) = smi_to_atoms_and_bonds(&smi).unwrap();
/// perceive_default_bonds(&atoms, &mut bonds);
/// assert_eq!(&bonds[0], &Bond::new(0, 1, BondType::Single));
/// ```
fn perceive_default_bonds(atoms: &Vec<Atom>, bonds: &mut Vec<Bond>) {
    let atom_idx_is_aromatic: Vec<bool> = atoms.iter().map(|atom| atom.aromatic).collect();
    for bond in bonds {
        if bond.bond_type == BondType::Default {
            if atom_idx_is_aromatic[bond.atom_idx_1] && atom_idx_is_aromatic[bond.atom_idx_2] {
                bond.bond_type = BondType::Aromatic;
            } else {
                bond.bond_type = BondType::Single;
            }
        }
    }
}

/// Gives the valence and indices of the neighbors of an atom specified by its atom_idx.
/// 
/// # Arguments
/// * `atom_idx`: The index of the atom in the molecule.
/// * `bonds`: A vector of all the bonds in a molecule.
/// 
/// # Examples
/// ```
/// let smi = "CC";
/// let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds(&smi).unwrap();
/// let (valence, neighbor_idxs) = get_valence_and_neighbor_idxs(0, &bonds);
/// assert_eq!(&valence, &1);
/// assert_eq!(&neighbor_idxs, &vec![1]);
/// ```
fn get_explicit_valence_and_neighbor_idxs(atom_idx: usize, bonds: &Vec<Bond>) -> (usize, Vec<usize>) {
    let mut valence = 0.0;
    let mut neighbor_idxs = vec![];

    for bond in bonds.iter() {
        if bond.atom_idx_1 == atom_idx || bond.atom_idx_2 == atom_idx {
            valence += bond.bond_type.to_float();
        }
        if bond.atom_idx_1 == atom_idx { neighbor_idxs.push(bond.atom_idx_2); }
        if bond.atom_idx_2 == atom_idx { neighbor_idxs.push(bond.atom_idx_1); }
    }

    (valence as usize, neighbor_idxs)
}

fn atom_is_double_bonded_to_o(atom_idx: usize, atoms: &Vec<Atom>, bonds: &Vec<Bond>) -> bool {
    for bond in bonds.iter() {
        if bond.bond_type == BondType::Double {
            if bond.atom_idx_1 == atom_idx {
                if atoms[bond.atom_idx_2].atomic_symbol == AtomicSymbol::O { return true; }
            } else if bond.atom_idx_2 == atom_idx {
                if atoms[bond.atom_idx_1].atomic_symbol == AtomicSymbol::O { return true; }
            }
        }
    }

    false
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_empty_smi() {
        let smi = "";
        let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
        assert_eq!(&atoms, &vec![]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn parse_linear_smi() {
        let smi = "C";
        let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            Atom::new(AtomicSymbol::C, 0, 0, false, 4, Chirality::Undefined, 0, 0, vec![], 0, "C"),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());
        let smi = "ClCCl";
        let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            Atom::new(AtomicSymbol::Cl, 0, 0, false, 0, Chirality::Undefined, 0, 0, vec![1], 0, "Cl"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 2, Chirality::Undefined, 0, 1, vec![0, 2], 0, "C"),
            Atom::new(AtomicSymbol::Cl, 0, 0, false, 0, Chirality::Undefined, 0, 2, vec![1], 0, "Cl"),
        ]);
        assert_eq!(&bonds, &vec![
            Bond::new(0, 1, BondType::Single),
            Bond::new(1, 2, BondType::Single),
        ]);
        assert_eq!(&ring_closures, &HashMap::new());
        let smi = "OCCOCCOCCO";
        let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            Atom::new(AtomicSymbol::O, 0, 0, false, 1, Chirality::Undefined, 0, 0, vec![1], 0, "O"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 2, Chirality::Undefined, 0, 1, vec![0, 2], 0, "C"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 2, Chirality::Undefined, 0, 2, vec![1, 3], 0, "C"),
            Atom::new(AtomicSymbol::O, 0, 0, false, 0, Chirality::Undefined, 0, 3, vec![2, 4], 0, "O"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 2, Chirality::Undefined, 0, 4, vec![3, 5], 0, "C"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 2, Chirality::Undefined, 0, 5, vec![4, 6], 0, "C"),
            Atom::new(AtomicSymbol::O, 0, 0, false, 0, Chirality::Undefined, 0, 6, vec![5, 7], 0, "O"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 2, Chirality::Undefined, 0, 7, vec![6, 8], 0, "C"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 2, Chirality::Undefined, 0, 8, vec![7, 9], 0, "C"),
            Atom::new(AtomicSymbol::O, 0, 0, false, 1, Chirality::Undefined, 0, 9, vec![8], 0, "O"),
        ]);
        assert_eq!(&bonds, &vec![
            Bond::new(0, 1, BondType::Single),
            Bond::new(1, 2, BondType::Single),
            Bond::new(2, 3, BondType::Single),
            Bond::new(3, 4, BondType::Single),
            Bond::new(4, 5, BondType::Single),
            Bond::new(5, 6, BondType::Single),
            Bond::new(6, 7, BondType::Single),
            Bond::new(7, 8, BondType::Single),
            Bond::new(8, 9, BondType::Single),
        ]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn parse_bond_smi() {
        let smi = "C-C";
    }

    #[test]
    fn parse_bracket_atom() {
        let smi = "Cl";
        let (atoms, bonds, ring_closures) = smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
        assert_eq!(&atoms, &vec![
            Atom::new(AtomicSymbol::Cl, 0, 0, false, 1, Chirality::Undefined, 0, 0, vec![], 0, "Cl"),
        ]);
        assert_eq!(&bonds, &vec![]);
        assert_eq!(&ring_closures, &HashMap::new());
    }

    #[test]
    fn test_get_explicit_valence_and_neighbor_idxs() {
        let bonds = vec![
            Bond::new(0, 1, BondType::Single),
            Bond::new(0, 2, BondType::Single),
            Bond::new(0, 3, BondType::Double),
        ];
        let (valence, neighbor_idxs) = get_explicit_valence_and_neighbor_idxs(0, &bonds);

        assert_eq!(valence, 4);
        assert_eq!(neighbor_idxs, vec![1, 2, 3]);
    }

    #[test]
    fn test_perceive_default_bonds() {
        let atoms = vec![
            Atom::new(AtomicSymbol::C, 0, 0, true, 0, Chirality::Undefined, 0, 0, vec![], 0, "C"),
            Atom::new(AtomicSymbol::C, 0, 0, true, 0, Chirality::Undefined, 0, 1, vec![], 0, "C"),
        ];
        let mut bonds = vec![
            Bond::new(0, 1, BondType::Default),
        ];
        perceive_default_bonds(&atoms, &mut bonds);
        assert_eq!(&bonds, &vec![Bond::new(0, 1, BondType::Aromatic),]);
        let atoms = vec![
            Atom::new(AtomicSymbol::C, 0, 0, true, 0, Chirality::Undefined, 0, 0, vec![], 0, "C"),
            Atom::new(AtomicSymbol::C, 0, 0, false, 0, Chirality::Undefined, 0, 1, vec![], 0, "C"),
        ];
        let mut bonds = vec![
            Bond::new(0, 1, BondType::Default),
        ];
        perceive_default_bonds(&atoms, &mut bonds);
        assert_eq!(&bonds, &vec![Bond::new(0, 1, BondType::Single),]);
    }
}