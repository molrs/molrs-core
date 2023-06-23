use std::{env, fs, process::Command};

use petgraph::dot::Dot;

use crate::molecule::Molecule;

pub fn sd_file() {}

pub fn pdb_file() {}

pub fn xyz_file() {}

pub fn dot_png_file(mol: &Molecule, path: &str) {
    let dot_str = format!("{:?}", Dot::with_config(&mol.graph, &[]));
    let tmp_dir = env::temp_dir();
    fs::write(format!("{}/graph.dot", tmp_dir.to_str().unwrap()), dot_str)
        .expect("unable to write graph.dot");
    Command::new("dot")
        .arg("-Tpng")
        .arg(format!("{}/graph.dot", tmp_dir.to_str().unwrap()))
        .arg("-o")
        .arg(path)
        .output()
        .unwrap();
}
