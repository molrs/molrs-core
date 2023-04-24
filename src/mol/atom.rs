/// TODO:
/// hybridization?

use crate::periodic_table::AtomicSymbol;
use super::Mol;
use super::Bond;
use super::bond::BondType;


#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Chirality {
    #[default]
    Undefined,
    Clockwise,
    CounterClockwise,
}

impl Chirality {
    pub fn new(chirality: &str) -> Result<Chirality, String> {
        if chirality == "" { Ok(Chirality::Undefined) }
        else if chirality == "@@" { Ok(Chirality::Clockwise) }
        else if chirality == "@" { Ok(Chirality::CounterClockwise) }
        else { return Err("invalid chirality".to_string()); }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Atom {
    pub atomic_symbol: AtomicSymbol,
    pub isotope: u16,
    pub charge: i8,
    pub aromatic: bool,
    pub num_imp_h: u8,
    pub chirality: Chirality,
    pub num_rad_electron: usize,
    pub idx: usize,
    pub neighbor_idxs: Vec<usize>,
    pub smallest_ring_size: usize,
}

impl Atom {
    pub fn new(
        atomic_symbol: AtomicSymbol,
        isotope: u16,
        charge: i8,
        aromatic: bool,
        num_imp_h: u8,
        chirality: Chirality,
        num_rad_electron: usize,
        idx: usize,
        neighbor_idxs: Vec<usize>,
        smallest_ring_size: usize,
    ) -> Atom {
        Atom {
            atomic_symbol: atomic_symbol,
            isotope: isotope,
            charge: charge,
            aromatic: aromatic,
            num_imp_h: num_imp_h,
            chirality: chirality,
            num_rad_electron: num_rad_electron,
            idx: idx,
            neighbor_idxs: neighbor_idxs,
            smallest_ring_size: smallest_ring_size,
        }
    }

    fn add_chirality(&self) {
        unimplemented!();
    }

    pub fn neighbors<'a>(&'a self, mol: &'a Mol) -> Vec<&'a Atom> {
        let mut neighbors: Vec<&Atom> = vec![];
        for i in &self.neighbor_idxs {
            neighbors.push(&mol.atoms[*i]);
        }

        neighbors
    }

    pub fn symbol(&self) -> String {
        unimplemented!();
    }

    pub fn atomic_number(&self) -> usize {
        unimplemented!();
    }

    pub fn atomic_mass(&self) -> f64 {
        unimplemented!();
    }

    pub fn bonds<'a>(&'a self, mol: &'a Mol) -> Vec<&'a Bond> {
        let mut bonds: Vec<&Bond> = vec![];
        for bond in &mol.bonds {
            if bond.atom_idx_1 == self.idx || bond.atom_idx_2 == self.idx { bonds.push(bond); }
        }

        bonds
    }

    pub fn explicit_valence<'a>(&'a self, mol:&'a Mol) -> usize {
        let mut valence = 0.0;
        for bond in &mol.bonds {
            valence += bond.bond_type.to_float();
        }

        valence as usize
    }

    pub fn total_valence<'a>(&'a self, mol: &'a Mol) -> usize {
        self.explicit_valence(&mol) + self.num_imp_h as usize
    }

    pub fn max_allowed_valence(&self) -> Option<usize> {
        if self.atomic_symbol == AtomicSymbol::B {
            if self.charge == 0 { Some(3) }
            else if self.charge == 1 { Some(2) }
            else if self.charge == -1 { Some(4) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::C ||
        self.atomic_symbol == AtomicSymbol::Si {
            if self.charge == 0 { Some(4) }
            else if self.charge == 1 || self.charge == -1 { Some(3) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::N {
            if self.charge == 0 { Some(3) }
            else if self.charge == 1 { Some(4) }
            else if self.charge == -1 { Some(2) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::O {
            if self.charge == 0 { Some(2) }
            else if self.charge == 1 { Some(3) }
            else if self.charge == -1 { Some(1) }
            else if self.charge == -2 { Some(0) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::F {
            if self.charge == 0 { Some(1) }
            else if self.charge == -1 { Some(0) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::P ||
        self.atomic_symbol == AtomicSymbol::As || self.atomic_symbol == AtomicSymbol::Sb ||
        self.atomic_symbol == AtomicSymbol::Bi {
            if self.charge == 0 { Some(5) }
            else if self.charge == 1 || self.charge == -1 { Some(4) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::S ||
        self.atomic_symbol == AtomicSymbol::Se || self.atomic_symbol == AtomicSymbol::Te {
            if self.charge == 0 { Some(6) }
            else if self.charge == 1 || self.charge == -1 { Some(5) } 
            else if self.charge == -2 { Some(4) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::Cl ||
        self.atomic_symbol == AtomicSymbol::Br || self.atomic_symbol == AtomicSymbol::I {
            if self.charge == 0 { Some(7) }
            else if self.charge == 1 || self.charge == -1 { Some(6) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::H ||
        self.atomic_symbol == AtomicSymbol::Li || self.atomic_symbol == AtomicSymbol::Na ||
        self.atomic_symbol == AtomicSymbol::K || self.atomic_symbol == AtomicSymbol::Rb ||
        self.atomic_symbol == AtomicSymbol::Cs || self.atomic_symbol == AtomicSymbol::Fr ||        
        self.atomic_symbol == AtomicSymbol::Ag {
            if self.charge == 0 { Some(1) }
            else if self.charge == 1 { Some(0) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::Be ||
        self.atomic_symbol == AtomicSymbol::Mg || self.atomic_symbol == AtomicSymbol::Ca ||
        self.atomic_symbol == AtomicSymbol::Sr || self.atomic_symbol == AtomicSymbol::Ba ||
        self.atomic_symbol == AtomicSymbol::Ra || self.atomic_symbol == AtomicSymbol::Zn {
            if self.charge == 0 { Some(2) }
            else if self.charge == 1 { Some(1) }
            else if self.charge == 2 { Some(0) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::Al ||
        self.atomic_symbol == AtomicSymbol::Fe {
            if self.charge == 0 { Some(3) }
            else if self.charge == 1 { Some(2) }
            else if self.charge == 2 { Some(1) }
            else if self.charge == 3 { Some(0) }
            else { None }
        } else if self.atomic_symbol == AtomicSymbol::Co {
            if self.charge == 0 { Some(5) }
            else if self.charge == 1 { Some(4) }
            else if self.charge == 2 { Some(3) }
            else if self.charge == 3 { Some(2) }
            else if self.charge == 4 { Some(1) }
            else if self.charge == 5 { Some(0) }
            else { None }
        } else {
            None
        }
    }

    pub fn explicit_degree(&self) -> usize {
        self.neighbor_idxs.len()
    }

    pub fn total_degree(&self) -> usize {
        self.neighbor_idxs.len() + self.num_imp_h as usize
    }

    pub fn total_num_h<'a>(&'a self, mol: &'a Mol) -> usize {
        let mut total_num_h = self.num_imp_h;
        for neighbor_idx in &self.neighbor_idxs {
            let neighbor = &mol.atoms[*neighbor_idx];
            if neighbor.atomic_symbol == AtomicSymbol::H { total_num_h += 1; }
        }

        total_num_h as usize
    }

    pub fn has_double_bond<'a>(&'a self, mol: &'a Mol) -> bool {
        self.bonds(mol)
            .iter()
            .map(|bond| &bond.bond_type)
            .filter(|bond_type| *bond_type == &BondType::Double)
            .count() > 0
    }

    pub fn is_conjugated<'a>(&'a self, mol: &'a Mol) -> bool {
        self.aromatic || self.has_double_bond(mol)
    }

    pub fn is_chiral_center<'a>(&'a self, mol: &'a Mol) -> bool {
        let mut num_double_bond = 0;
        for bond in self.bonds(mol) {
            if bond.bond_type == BondType::Double { num_double_bond += 1; }
        }
        if num_double_bond == 2 { return self.is_allene_chiral_center(mol) }
        else { return self.is_tet_chiral_center(mol) }
    }

    fn is_allene_chiral_center<'a>(&'a self, mol: &'a Mol) -> bool {
        true
    }

    fn is_tet_chiral_center<'a>(&'a self, mol: &'a Mol) -> bool {
        if self.total_num_h(mol) > 1 || self.num_imp_h as usize + self.neighbor_idxs.len() < 4 { return false; }

        let mut neighbor_paths: Vec<Box<Vec<usize>>> = vec![];
        for start_idx in &self.neighbor_idxs {
            let start_idx = *start_idx;
            neighbor_paths.push(Box::new(vec![self.idx, start_idx]));
        }

        let mut new_neighbor_paths: Vec<Box<Vec<usize>>> = vec![];
        for rank in 0..5 as usize {
            // instead of building a tree, try to chase down similarities?
            // i think we'll just have to build a tree though :(
        }

        true
    }
}