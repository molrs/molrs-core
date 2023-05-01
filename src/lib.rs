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
    // use atom::Atom;
    // use bond::Bond;
    use petgraph::dot::{Config, Dot};
    use write;
    // use petgraph::graph::EdgeReference;
    // use petgraph::visit::EdgeRef;
    use petgraph::{Graph, Undirected};
    use std::fs;

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
        // let smi = "[CH3]";
        // let smi = "C[CH2-]C";
        // let smi = "CCC";
        // let smi = "C1CC1";
        // let smi = "c1occc1C";
        let smi = "c12ncccc1[nH]cc2";
        // let smi = "c1%02ncccc1[nH]cc%02";
        // let smi = "C1OC=CC=1";
        // let smi = "C1C[CH2]1";
        // let smi = "C%01CC%01";
        // let smi = "C1CC=1";
        // let smi = "C1CC=1-C";
        // let smi = "C%01CC=%01";
        let mol = from::smiles(smi).unwrap();
        dbg!(&mol);
        write::dot_png_file(&mol, "graph.png");
    }

    #[test]
    fn test_atom() {
        // let mut graph: Graph<Atom, Bond, Undirected> = Graph::new_undirected();
        // graph.add_node(Atom{idx: graph.node_count(), atomic_number: 6});
        // graph.add_node(Atom{idx: graph.node_count(), atomic_number: 6});
        // graph.add_edge(0.into(), 1.into(), Bond{atom_idx_1: 0, atom_idx_2: 1, bond_type: bond::BondType::Default});
        // dbg!(&graph);
        // dbg!(graph.edges(0.into()).collect::<Vec<EdgeReference<Bond>>>()[0].source());
    }

    #[test]
    fn test_graph() {
        let mut graph: Graph<u8, u8, Undirected> = Graph::new_undirected();
        graph.add_node(1);
        graph.add_node(1);
        graph.add_edge(0.into(), 1.into(), 1);
        let dot_str = format!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
        fs::write("graph.dot", dot_str).expect("unable to write graph.dot");
    }
}
