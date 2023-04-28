#![feature(test)]

pub mod mol;
pub mod sub;
pub mod read;
pub mod periodic_table;

use mol::Mol;
use sub::Sub;

#[cfg(test)]
mod tests {
    extern crate test;

    use std::fs::File;
    use std::io::{self, BufRead};
    use std::path::Path;
    use std::time::{Duration, Instant};
    use test::{Bencher, black_box};
    use super::*;


    // fn smi_to_atoms_and_bonds_no_return(smi: &str) {
    //     let (atoms, bonds, ring_closures) = read::smiles::smi_to_atoms_and_bonds_and_ring_closures(&smi).unwrap();
    // }

    // fn mol_from_smiles_no_return(smi: &str) {
    //     let mol = Mol::from_smiles(&smi).unwrap();
    // }

    // #[test]
    // fn check_smi() {
    //     // let smi = "N[Cu](N)(N)[OH-]";
    //     // let smi = "N~[Cu](~N)(~[OH-])~[OH-]";
    //     // let smi = "CS(=NC)(=O)C";
    //     // let smi = "o1cccc1";
    //     let smi = "CC(C)(C)(C)C";
    //     let mol = Mol::from_smiles(&smi);
    //     dbg!(&mol);
    // }

    // #[bench]
    // fn benchmark(b: &mut Bencher) {
    //     // let smi = "OCC[OH]";
    //     let smi = "CC(C)(C#CC1=NC(=C(C=C1)C2=C3C(=C(C=C2)Cl)C(=NN3CC(F)(F)F)NS(=O)(=O)C)C(CC4=CC(=CC(=C4)F)F)NC(=O)CN5C6=C(C7CC7C6(F)F)C(=N5)C(F)(F)F)S(=O)(=O)C";
    //     b.iter(|| {
    //         black_box(mol_from_smiles_no_return(&smi));
    //         // black_box(smi_to_atoms_and_bonds_no_return(&smi));
    //     });
    // }

    // fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    // where P: AsRef<Path>, {
    //     let file = File::open(filename)?;
    //     Ok(io::BufReader::new(file).lines())
    // }

    // base: 2.17 s
    // perceive rings: 9.55 s
    // perceive rings + kekulize: 10.48 s
    // #[test]
    // fn read_chembl_mols() {
    //     let mut smiles = vec![];
    //     if let Ok(lines) = read_lines("/Users/ozone/Downloads/chembl.smi") {
    //         for line in lines {
    //             if let Ok(ip) = line {
    //                 let smi = ip.split(" ").nth(0).unwrap();
    //                 smiles.push(smi.to_owned());
    //             }
    //         }
    //     }
        
    //     let start = Instant::now();

    //     for smi in smiles {
    //         match Mol::from_smiles(&smi) {
    //             Ok(_) => (),
    //             Err(problem) => println!("{}", &problem),
    //             // Err(_) => println!("{}", &smi),
    //         }
    //         // let result = Mol::from_smiles(&smi);
    //         // let result = read::smiles::smi_to_atoms_and_bonds_and_ring_closures(&smi);
    //     }

    //     let duration = start.elapsed();
    //     println!("Time elapsed is: {:?}", duration);
    // }
}
