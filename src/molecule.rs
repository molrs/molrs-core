use std::str::FromStr;

use crate::{
    atom::Atom,
    bond::Bond,
    smiles::{SmilesParseError, SmilesParser},
};

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
}

impl FromStr for Molecule {
    type Err = SmilesParseError;

    fn from_str(smi: &str) -> Result<Self, Self::Err> {
        let smiles_parser: SmilesParser = smi.parse()?;
        let atoms = smiles_parser.atoms()?;
        let bonds = smiles_parser.bonds()?;

        // perceive rings
        // perceive default bonds
        // kekulize
        // perceive implicit hydrogens
        // delocalize

        Ok(Molecule { atoms, bonds })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_molecule_from_str() {
        let smi = "C2CC2";
        dbg!(Molecule::from_str(smi).unwrap());
    }
}
