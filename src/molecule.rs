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
