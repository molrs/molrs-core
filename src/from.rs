use std::collections::HashMap;

use crate::{atom::Atom, bond::BondType, molecule::Molecule};

fn parse_smiles(smi: &str) -> (Vec<String>, Vec<String>, Vec<String>) {
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

    (atom_strs, bond_strs, ring_closure_strs)
}

fn add_atoms_to_molecule(molecule: &mut Molecule, atom_strs: &[String]) {
    atom_strs.iter().for_each(|atom_str| {
        molecule
            .graph
            .add_node(Atom::from_str(atom_str, molecule.graph.node_count()).unwrap());
    });
}

fn add_non_ring_closure_bonds_to_molecule(
    molecule: &mut Molecule,
    bond_strs: &[String],
) -> Result<(), ()> {
    let mut source_atom_index_stack = vec![0];
    for (i, bond_str) in bond_strs.iter().enumerate() {
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
                return Err(());
            }
        }
        if bond_char != '.' {
            molecule.graph.add_edge(
                source_atom_index_stack.pop().unwrap().into(),
                target_atom_index.into(),
                match BondType::from_char(bond_char) {
                    Ok(bond_type) => bond_type,
                    Err(_) => {
                        return Err(());
                    }
                },
            );
            source_atom_index_stack.push(target_atom_index);
        }
    }

    Ok(())
}

#[derive(Debug)]
enum RingClosure {
    OneAtom(usize),
    Used,
}

fn add_ring_closure_bonds_to_molecule(
    molecule: &mut Molecule,
    ring_closure_strs: &[String],
) -> Result<(), ()> {
    let mut ring_closures = HashMap::new();
    for (i, ring_closure_str) in ring_closure_strs.iter().enumerate() {
        if ring_closure_str.is_empty() {
            continue;
        }
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
                let last_element = ring_indices.last_mut().unwrap();
                *last_element = *last_element * 10 + c.to_digit(10).unwrap();
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
                return Err(());
            }
        }

        for ring_index in ring_indices {
            match ring_closures.get(&ring_index) {
                Some(ring_closure) => match ring_closure {
                    RingClosure::OneAtom(ring_closure_atom_index) => {
                        molecule.graph.add_edge(
                            (*ring_closure_atom_index).into(),
                            source_atom_index.into(),
                            BondType::from_char(bond_char).unwrap(),
                        );
                        ring_closures.insert(ring_index, RingClosure::Used);
                    }
                    RingClosure::Used => {
                        return Err(());
                    }
                },
                None => {
                    ring_closures.insert(ring_index, RingClosure::OneAtom(source_atom_index));
                }
            }
        }
    }

    Ok(())
}

pub fn smiles(smi: &str) -> Result<Molecule, String> {
    if smi.chars().next().is_none() {
        return Err("SMILES string is empty".to_owned());
    }

    let mut molecule = Molecule::default();

    let (atom_strs, bond_strs, ring_closure_strs) = parse_smiles(smi);

    add_atoms_to_molecule(&mut molecule, &atom_strs);
    match add_non_ring_closure_bonds_to_molecule(&mut molecule, &bond_strs) {
        Ok(_) => (),
        Err(_) => return Err(format!("error in smi {}", &smi)),
    };
    match add_ring_closure_bonds_to_molecule(&mut molecule, &ring_closure_strs) {
        Ok(_) => (),
        Err(_) => return Err(format!("error in smi {}", &smi)),
    };

    molecule.perceive_rings();
    molecule.perceive_default_bonds();
    molecule = match molecule.kekulized() {
        Ok(mol) => mol,
        Err(err) => {
            println!("{}", &err);
            return Err(format!("kekulization failed for smi {}", &smi));
        }
    };
    match molecule.perceive_implicit_hydrogens() {
        Ok(_) => (),
        Err(_) => return Err(format!("valence problems in smi {}", &smi)),
    };
    // molecule.perceive_radicals();
    // molecule = match molecule.aromatized() {
    //     Ok(mol) => mol,
    //     Err(_) => return Err(format!("aromatizaton failed for smi {}", &smi)),
    // };

    Ok(molecule)
}

pub fn smarts() {}

pub fn mdl() {}

pub fn pdb() {}

pub fn xyz() {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_smiles() {
        let smiles = [
            "C",
            "CC",
            "C=C",
            "CC=O",
            "C[H]",
            "[18O]C",
            "C1CC=1",
            "C%10C[N-]%10",
            "C(F)C",
            "C(C)(C)C",
            "C(C(C)C)C",
        ];
        let vec_tuples = [
            (vec!["C"], vec![], vec![""]),
            (vec!["C", "C"], vec![""], vec!["", ""]),
            (vec!["C", "C"], vec!["="], vec!["", ""]),
            (vec!["C", "C", "O"], vec!["", "="], vec!["", "", ""]),
            (vec!["C", "[H]"], vec![""], vec!["", ""]),
            (vec!["[18O]", "C"], vec![""], vec!["", ""]),
            (vec!["C", "C", "C"], vec!["", ""], vec!["1", "", "=1"]),
            (vec!["C", "C", "[N-]"], vec!["", ""], vec!["%10", "", "%10"]),
            (vec!["C", "F", "C"], vec!["(", ")"], vec!["", "", ""]),
            (
                vec!["C", "C", "C", "C"],
                vec!["(", ")(", ")"],
                vec!["", "", "", ""],
            ),
            (
                vec!["C", "C", "C", "C", "C"],
                vec!["(", "(", ")", ")"],
                vec!["", "", "", "", ""],
            ),
        ];

        for (smi, vec_tuple) in smiles.iter().zip(vec_tuples) {
            let atom_strs = vec_tuple.0.iter().map(|x| String::from(*x)).collect();
            let bond_strs = vec_tuple.1.iter().map(|x| String::from(*x)).collect();
            let ring_closure_strs = vec_tuple.2.iter().map(|x| String::from(*x)).collect();
            assert_eq!(parse_smiles(smi), (atom_strs, bond_strs, ring_closure_strs));
        }
    }
}
