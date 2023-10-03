// use svg::node::element::path::Data;
// use svg::node::element::Path;
// use svg::Document;

// use crate::molecule::Molecule;

// pub fn svg_str(mol: &Molecule) -> String {
//     let data = Data::new()
//         .line_to((0, 10));

//     let path = Path::new()
//         .set("fill", "none")
//         .set("stroke", "black")
//         .set("stroke-width", 3)
//         .set("d", data);

//     let document = Document::new().set("viewBox", (0, 0, 70, 70)).add(path);

//     document.to_string()
// }

pub fn smiles() {}

pub fn canonical_smiles() {}

pub fn smarts() {}

pub fn mdl() {}

pub fn pdb() {}

pub fn xyz() {}

#[cfg(test)]
mod tests {
    // use super::*;
    // use std::fs;

    // #[test]
    // fn test_to_image() {
    //     let mol = crate::from::smiles("CC").unwrap();
    //     fs::write("mol.svg", svg_str(&mol)).expect("Unable to write file");
    // }
}
