/// TODO:
/// check for mol correctness? should be correct after read smiles but is needed when reading mol blocks
/// write Vec<Mol> to sdf and read Vec<Mol> from sdf
/// writing and reading user-defined properties?
/// running reactions
/// morgan bits and fingerprint
/// clustering
/// descriptors
/// fragmenting methods, ie fragment to obtain rings, fragment on synthetically accessible bonds

pub mod atom;
pub mod bond;
pub mod ring;
pub mod conf;

use std::collections::HashMap;

use atom::Atom;
use bond::Bond;
use ring::Ring;
use conf::Conf;
use self::{atom::Chirality, bond::BondType};

use super::sub::Sub;
use crate::{read::smiles, periodic_table::AtomicSymbol};

#[derive(Debug, Clone)]
pub struct Mol {
    name: String,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    rings: Vec<Ring>,
    confs: Vec<Conf>,

    // maybe don't need to be stored?
    // num_heavy_atoms: usize,
}

impl Mol {
    pub fn new() -> Mol {
        Mol {
            name: "".to_owned(),
            atoms: vec![],
            bonds: vec![],
            rings: vec![],
            confs: vec![],
        }
    }

    pub fn from_mol(mol: &Mol) -> Mol {
        // this function may not be necessary?
        Mol {
            name: mol.name.clone(),
            atoms: mol.atoms.clone(),
            bonds: mol.bonds.clone(),
            rings: mol.rings.clone(),
            confs: mol.confs.clone(),
        }
    }

    pub fn from_smiles(smi: &str) -> Result<Mol, String> {
        let (atoms, bonds, ring_closures) = match smiles::smi_to_atoms_and_bonds_and_ring_closures(&smi) {
            Ok((atoms, bonds, ring_closures)) => (atoms, bonds, ring_closures),
            Err(string) => return Err(string),
        };

        let mut mol = Mol {
            name: smi.to_owned(),
            atoms: atoms,
            bonds: bonds,
            rings: vec![],
            confs: vec![],
        };

        // convert default bonds to single or aromatic
        // convert L-type metal - non-metal bonds to dative bonds
        // break X-type metal - non-metal bonds and charge separate
        // explicit valence check
        mol.perceive_rings_from_ring_closures(&ring_closures);
        // mol.perceive_unnecessary_stereochem();
        // mol.add_2d_coordinates();
        mol = match mol.to_kekulized() {
            Ok(mol) => mol,
            Err(string) => return Err(string),
        };
        // set implicit hydrogens
        // canonicalize atom order

        Ok(mol)
    }

    pub fn from_mol_block() {
        unimplemented!();
    }

    fn perceive_rings_from_ring_closures(&mut self, ring_links: &HashMap<usize, (usize, usize)>) {
        for ring_link in ring_links.values() {
            let atom_idx_1 = ring_link.0;
            let atom_idx_2 = ring_link.1 as usize;
            
            let mut ring_paths: Vec<Box<Vec<usize>>> = vec![];
            for start_idx in &self.atoms[atom_idx_1].neighbor_idxs {
                let start_idx = *start_idx;
                if start_idx != atom_idx_2 {
                    ring_paths.push(Box::new(vec![atom_idx_1, start_idx]));
                }
            }

            loop {
                let mut new_ring_paths: Vec<Box<Vec<usize>>> = vec![];
                let mut reached_atom_2 = false;
                for path in ring_paths.iter() {
                    let path_head = *path.iter().rev().next().unwrap();
                    if path_head == atom_idx_2 {
                        new_ring_paths.push(path.clone());
                        reached_atom_2 = true;
                    } else {
                        for neighbor in &self.atoms[path_head].neighbor_idxs {
                            if !path.contains(&neighbor) {
                                let mut new_path = path.clone();
                                new_path.push(*neighbor);
                                new_ring_paths.push(new_path);
                            }
                        }
                    }
                }
                ring_paths = new_ring_paths;
                if reached_atom_2 {
                    // pop ring paths that don't end with atom_idx_2
                    let mut new_ring_paths = vec![];
                    for _ in 0..ring_paths.len() {
                        let ring_path = ring_paths.pop().unwrap();
                        if *ring_path.last().unwrap() == atom_idx_2 {
                            new_ring_paths.push(ring_path)
                        }
                    }
                    ring_paths = new_ring_paths;
                    break;
                } else if ring_paths.len() == 0 {
                    break;
                }
            }

            for ring_path in ring_paths.iter_mut() {
                let len_ring_path = ring_path.len();
                for atom_idx in ring_path.iter() {
                    let smallest_ring_size = self.atoms[*atom_idx].smallest_ring_size;
                    if smallest_ring_size == 0 || len_ring_path < smallest_ring_size {
                        self.atoms[*atom_idx].smallest_ring_size = len_ring_path;
                    }
                }
                let new_ring = Ring::new(&ring_path);
                let mut new_ring_is_in_rings = false;
                for ring in &self.rings {
                    if new_ring.atom_idxs.len() != ring.atom_idxs.len() { continue; }
                    let mut num_atom_idx_matches = 0;
                    for atom_idx_i in &new_ring.atom_idxs {
                        for atom_idx_j in &ring.atom_idxs {
                            if atom_idx_i == atom_idx_j {
                                num_atom_idx_matches += 1;
                            }
                        }
                    }
                    if num_atom_idx_matches == new_ring.atom_idxs.len() {
                        new_ring_is_in_rings = true;
                    }
                }
                if !new_ring_is_in_rings { self.rings.push(new_ring) }

                // if !self.rings.contains(&ring) {
                //     self.rings.push(Ring::new(&ring_path));
                // }
            }

            self.rings.sort_by_key(|k| k.atom_idxs.len());
        }
    }

    fn perceive_unnecessary_stereochem(&mut self) {
        unimplemented!();

        let mut stereochems_are_unnecessary: Vec<bool> = vec![];
        for atom in &self.atoms {
            if [Chirality::Undefined, Chirality::Clockwise, Chirality::CounterClockwise].contains(&atom.chirality) {
                if !atom.is_chiral_center(self) {
                    stereochems_are_unnecessary.push(true);
                    continue;
                } else { stereochems_are_unnecessary.push(false); }
            }
        }

        for (atom, stereochem_is_unnecessary) in self.atoms.iter_mut().zip(stereochems_are_unnecessary.iter()) {
            if *stereochem_is_unnecessary { atom.chirality = Chirality::Undefined; }
        }
    }

    fn add_2d_coordinates(&mut self) {
        // is automatically called when creating new mol
        // adds 2d coordinates in place
        unimplemented!();
    }

    pub fn to_smiles(&self) {
        unimplemented!();
    }

    pub fn to_mol_block(&self) {
        unimplemented!();
    }

    pub fn to_kekulized(&self) -> Result<Mol, String> {
        let mut mol = Mol {
            name: self.name.clone(),
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone().into_iter().filter(|bond| bond.bond_type != BondType::Aromatic).collect(),
            rings: self.rings.clone(),
            confs: self.confs.clone(),
        };

        let mut conjugated_rings = vec![];
        for ring in &self.rings {
            let mut ring_is_conjugated = true;
            for atom in ring.atom_idxs.iter().rev().map(|atom_idx| &self.atoms[*atom_idx]) {
                if !atom.is_conjugated(&mol) {
                    ring_is_conjugated = false;
                    break;
                }
            }
            if ring_is_conjugated { conjugated_rings.push(ring); }
        }

        let mut fused_rings_ring_idxs: Vec<Box<Vec<usize>>> = vec![];
        for (i, ring_i) in conjugated_rings.iter().enumerate() {
            for (j, ring_j) in conjugated_rings[(i+1)..].iter().enumerate() {
                let j = j + i + 1;
                if ring_i.atom_idxs.iter().filter(|atom_idx| ring_j.atom_idxs.contains(&atom_idx)).count() > 0 {
                    let mut part_of_existing_ring_system = false;
                    for fused_ring in fused_rings_ring_idxs.iter_mut() {
                        if fused_ring.contains(&i) {
                            fused_ring.push(j);
                            part_of_existing_ring_system = true;
                        } else if fused_ring.contains(&j) {
                            fused_ring.push(i);
                            part_of_existing_ring_system = true;
                        }
                    }
                    if !part_of_existing_ring_system {
                        fused_rings_ring_idxs.push(Box::new(vec![i, j]));
                    }
                }
            }
        }

        for i in 0..conjugated_rings.len() {
            let mut ring_is_fused = false;
            for fused_ring_ring_idxs in &fused_rings_ring_idxs {
                if fused_ring_ring_idxs.contains(&i) {
                    ring_is_fused = true;
                }
            }
            if !ring_is_fused {
                fused_rings_ring_idxs.push(Box::new(vec![i]));
            }
        }

        let mut fused_conjugated_rings = vec![];
        for fused_ring_ring_idxs in &fused_rings_ring_idxs {
            let mut atom_idxs = vec![];
            for ring in fused_ring_ring_idxs.iter().map(|ring_idx| conjugated_rings[*ring_idx]) {
                for atom_idx in &ring.atom_idxs {
                    if !atom_idxs.contains(atom_idx) {
                        atom_idxs.push(*atom_idx);
                    }
                }
            }
            fused_conjugated_rings.push(Ring::new(&atom_idxs));
        }

        for ring in &fused_conjugated_rings {
            let mut atoms_needing_double_bonds = vec![];
            for idx in &ring.atom_idxs {
                let atom = &self.atoms[*idx];
                if atom.aromatic &&
                atom.total_degree() < atom.max_allowed_valence().unwrap() &&
                atom.bonds(self).iter().map(|bond| &bond.bond_type).filter(|bond_type| *bond_type == &BondType::Double).count() == 0 {
                    atoms_needing_double_bonds.push(idx);
                } else {
                    if atom.aromatic {
                        let mut atom = mol.atoms.remove(*idx);
                        atom.aromatic = false;
                        mol.atoms.insert(*idx, atom);
                    }
                }
            }

            let mut continuous_chains: Vec<Box<Vec<usize>>> = vec![];
            while atoms_needing_double_bonds.len() != 0 {
                let mut continuous_chain = Box::new(vec![*atoms_needing_double_bonds.pop().unwrap()]);
                let mut neighbor_idxs = vec![];
                loop {
                    let len_continous_chain = continuous_chain.len();
                    for atom_idx in continuous_chain.iter() {
                        for neighbor_idx in &self.atoms[*atom_idx].neighbor_idxs {
                            neighbor_idxs.push(*neighbor_idx);
                        }
                    }
                    for neighbor_idx in &neighbor_idxs {
                        if !continuous_chain.contains(neighbor_idx) && atoms_needing_double_bonds.contains(&neighbor_idx) {
                            let index = atoms_needing_double_bonds.iter().position(|ele| *ele == neighbor_idx).unwrap();
                            continuous_chain.push(*atoms_needing_double_bonds.remove(index));
                        }
                    }
                    if len_continous_chain == continuous_chain.len() {
                        break;
                    }
                }
                continuous_chains.push(continuous_chain);
            }

            for continous_chain in &continuous_chains {
                if continous_chain.len() % 2 != 0 { return Err(format!("atom idxs {:?} could not be kekulized in smi {}", &continous_chain, &"")); }
                let mut chain_has_terminus = false;
                let mut terminus_idx = 0;
                for atom_idx in continous_chain.iter() {
                    if self.atoms[*atom_idx].neighbor_idxs.iter().filter(|neighbor_idx| continous_chain.contains(*neighbor_idx)).count() < 2 {
                        chain_has_terminus = true;
                        terminus_idx = *atom_idx;
                    }
                }
                let mut current_atoms;
                if chain_has_terminus {
                    current_atoms = vec![*continous_chain.iter().nth(continous_chain.iter().position(|ele| *ele == terminus_idx).unwrap()).unwrap()]
                } else {
                    current_atoms = vec![*continous_chain.first().unwrap()];
                }
                let mut traversed_atom_idxs = vec![];
                while traversed_atom_idxs.len() < continous_chain.len() {
                    let mut neighbor_idxs = vec![];
                    for atom_idx in &current_atoms {
                        traversed_atom_idxs.push(*atom_idx);
                        for neighbor_idx in &self.atoms[*atom_idx].neighbor_idxs {
                            if !traversed_atom_idxs.contains(&neighbor_idx) && continous_chain.contains(&neighbor_idx) {
                                neighbor_idxs.push(*neighbor_idx);
                                if mol.atoms[*atom_idx].has_double_bond(&mol) {
                                    mol.bonds.push(Bond::new(*atom_idx, *neighbor_idx, BondType::Single));
                                } else {
                                    mol.bonds.push(Bond::new(*atom_idx, *neighbor_idx, BondType::Double));
                                }
                            }
                        }
                        let mut atom = mol.atoms.remove(*atom_idx);
                        atom.aromatic = false;
                        mol.atoms.insert(*atom_idx, atom);
                    }
                    current_atoms = neighbor_idxs;
                }
            }
        }

        for atom in &mol.atoms {
            if atom.aromatic { return Err(format!("aromatic atom remaining in kekulized mol")) }
        }

        Ok(mol)
    }

    pub fn to_aromatized(&self) -> Result<Mol, String> {
        unimplemented!();
    }

    pub fn to_murcko_scaffold(&self) -> Mol {
        unimplemented!();
    }

    pub fn with_2d_coordinates(&self, template: Sub, partial: bool) {
        // allow generation of 2d coordinates with matching template
        unimplemented!();
    }

    pub fn to_svg(&self) {
        unimplemented!();
    }

    pub fn with_explicit_hydrogens(&self) -> Mol {
        let mut mol = self.clone();

        for (i, atom) in self.atoms.iter().enumerate() {
            for _ in 0..atom.num_imp_h {
                let num_atoms = mol.atoms.len();
                mol.atoms[i].num_imp_h -= 1;
                mol.atoms[i].neighbor_idxs.push(num_atoms);
                mol.bonds.push(Bond::new(
                    i,
                    num_atoms,
                    BondType::Single,
                ));
                mol.atoms.push(Atom::new(
                    AtomicSymbol::H,
                    0,
                    0,
                    false,
                    0,
                    Chirality::Undefined,
                    0,
                    num_atoms,
                    vec![i],
                    0,
                    "[H]",
                ));
            }
        }

        // set coordinates of new Hs

        mol
    }

    pub fn without_explicit_hydrogens(&self) -> Mol {
        let mut mol = self.clone();

        let mut atom_idxs_to_exclude: Vec<usize> = vec![];
        for (i, atom) in self.atoms.iter().enumerate() {
            if atom.atomic_symbol == AtomicSymbol::H {
                // TODO: this part could be cleaner
                let neighbor_idx = atom.neighbor_idxs[0];
                let mut neighbor_idxs: Vec<usize> = vec![];
                for j in mol.atoms[neighbor_idx].neighbor_idxs.iter() {
                    let j = *j;
                    if j != i { neighbor_idxs.push(j); }
                }
                mol.atoms[neighbor_idx].neighbor_idxs = neighbor_idxs;
                mol.atoms[neighbor_idx].num_imp_h += 1;
                atom_idxs_to_exclude.push(i)
            }
        }
        let mut bond_idxs_to_exclude: Vec<usize> = vec![];
        for (i, bond) in self.bonds.iter().enumerate() {
            if atom_idxs_to_exclude.contains(&bond.atom_idx_1) || atom_idxs_to_exclude.contains(&bond.atom_idx_2) {
                bond_idxs_to_exclude.push(i);
            }
        }

        let new_atoms = mol.atoms.iter().enumerate().filter(|i_and_atom| !atom_idxs_to_exclude.contains(&i_and_atom.0)).map(|i_and_atom| i_and_atom.1.clone()).collect();
        let new_bonds = mol.bonds.iter().enumerate().filter(|i_and_bond| !bond_idxs_to_exclude.contains(&i_and_bond.0)).map(|i_and_bond| i_and_bond.1.clone()).collect();
        // edit conf to remove hydrogens
        mol.atoms = new_atoms;
        mol.bonds = new_bonds;

        mol
    }

    pub fn atom_with_idx(&self, idx: usize) -> Result<&Atom, String> {
        match self.atoms.iter().nth(idx) {
            Some(atom) => Ok(atom),
            None => Err(format!("idx = {} is greater than self.atoms.len() = {}", idx, self.atoms.len())),
        }
    }

    pub fn bond_with_idx(&self, idx: usize) -> Result<&Bond, String> {
        match self.bonds.iter().nth(idx) {
            Some(bond) => Ok(bond),
            None => Err(format!("idx = {} is greater than self.bonds.len() = {}", idx, self.atoms.len())),
        }
    }

    pub fn bond_between_atoms(&self, atom_idx_1: usize, atom_idx_2: usize) -> Option<&Bond> {
        if atom_idx_1 == atom_idx_2 { return None; }

        let i: usize;
        let j: usize;
        if atom_idx_1 < atom_idx_2 {
            i = atom_idx_1;
            j = atom_idx_2;
        } else {
            i = atom_idx_2;
            j = atom_idx_1;
        }

        for bond in &self.bonds {
            if bond.atom_idx_1 == i && bond.atom_idx_2 == j { return Some(&bond); }
        }
        None
    }

    pub fn without_conformers(&self) -> Mol {
        Mol {
            name: self.name.clone(),
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
            rings: self.rings.clone(),
            confs: vec![],
        }
    }
    
    pub fn maximum_common_substructure(&self, other: &Mol) -> Mol {
        unimplemented!();
    }

    pub fn common_markush(&self, other: &Mol) -> Sub {
        unimplemented!();
    }

    pub fn match_mol(&self, other: &Mol) -> Vec<Vec<usize>> {
        unimplemented!();
    }

    pub fn match_sub(&self, other: &Sub) -> Vec<Vec<usize>> {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn playground() {
        // let smi = "C12CCC3C(CC3)CCC1C2";
        // let smi = "C[C@H](C)C([H])[H]";
        // let smi = "CS(=O)(=O)";
        // let smi = "c1ccccc1";
        // let smi = "C1=Ccccc1";
        // let smi = "c1ccc2ccccc2c1";
        // let smi = "c1ccc(Cc2ccc3ccccc3c2)cc1";
        // let smi = "c1ccc2cc(Cc3ccc4ccccc4c3)ccc2c1";
        let smi = "o1ccc(cccc2)c12";
        // let smi = "c1cccc(CCC2)c12";
        // let smi = "c1cc2c(ccc3ccoc32)[nH]1";
        // let smi = "O=c1ccn2ccccn12";
        let mol = Mol::from_smiles(&smi).unwrap();
        // let mol = mol.with_explicit_hydrogens();
        // let mol = mol.without_explicit_hydrogens();
        // dbg!(&mol);
    }
}