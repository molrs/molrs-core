#![feature(test)]

pub mod atom;
pub mod bond;
pub mod extended_atom;
pub mod from;
pub mod molecule;
pub mod periodic_table;
pub mod reaction;
pub mod read;
pub mod standardize;
pub mod substructure;
pub mod to;
pub mod utils;
pub mod write;

pub fn visualize() {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_smiles() {
        // let smi = "C";
        // let smi = "CN";
        // let smi = "C-C";
        // let smi = "C=C";
        // let smi = "CCNCCOCC";
        // let smi = "CC(C)CC";
        // let smi = "C(C(C)C)C";
        // let smi = "C(C)(C(C)C)C";
        // let smi = "C(F)(F)(F)F";
        // let smi = "CS=O";
        // let smi = "CCl";
        // let smi = "ClC";
        let smi = "[CH3-]";
        // let smi = "C[CH2-]C";
        // let smi = "CCC";
        // let smi = "C1CC1";
        // let smi = "c1occc1C";
        // let smi = "c1(=O)cc[nH]cc1";
        // let smi = "c1ccccc1";
        // let smi = "c1=ccccc1";
        // let smi = "c1c=cccc1";
        // let smi = "c12ncccc1[nH]cc2";
        // let smi = "c3(cco4)c4CCC3";
        // let smi = "c1%02ncccc1[nH]cc%02";
        // let smi = "C1OC=CC=1";
        // let smi = "C1C[CH2]1";
        // let smi = "C%01CC%01";
        // let smi = "C1CC=1";
        // let smi = "C1CC=1-C";
        // let smi = "C%01CC=%01";
        let mol = from::smiles(smi).unwrap();
        dbg!(&mol);
        // write::dot_png_file(&mol, "graph.png");
    }
}
