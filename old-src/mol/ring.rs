#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Ring {
    pub atom_idxs: Vec<usize>,
}

impl Ring {
    pub fn new(atom_idxs: &Vec<usize>) -> Ring {
        Ring {
            atom_idxs: atom_idxs.clone()
        }
    }
}