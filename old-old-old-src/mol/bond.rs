#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum BondType {
    #[default]
    Default,
    Single, Up, Down,
    Aromatic,
    Double,
    Triple,
    Quadruple,
    Dative,
}

impl BondType {
    pub fn from_char(bond_type: char) -> Option<BondType> {
        if bond_type == ' ' { Some(BondType::Default ) }
        else if bond_type == '-' { Some(BondType::Single) }
        else if bond_type == '/' { Some(BondType::Up) }
        else if bond_type == '\\' { Some(BondType::Down) }
        else if bond_type == ':' { Some(BondType::Aromatic) }
        else if bond_type == '=' { Some(BondType::Double) }
        else if bond_type == '#' { Some(BondType::Triple) }
        else if bond_type == '$' { Some(BondType::Quadruple) }
        else { None }
    }

    pub fn to_float(&self) -> f64 {
        match self {
            BondType::Default => 1.0,
            BondType::Single => 1.0,
            BondType::Up => 1.0,
            BondType::Down => 1.0,
            BondType::Aromatic => 1.5,
            BondType::Double => 2.0,
            BondType::Triple => 3.0,
            BondType::Quadruple => 4.0,
            BondType::Dative => 0.0,
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Bond {
    pub atom_idx_1: usize,
    pub atom_idx_2: usize,
    pub bond_type: BondType,
    pub smallest_ring_size: usize,
}

impl Bond {
    pub fn new(atom_idx_1: usize, atom_idx_2: usize, bond_type: BondType, smallest_ring_size: usize) -> Bond {
        Bond {
            atom_idx_1: atom_idx_1,
            atom_idx_2: atom_idx_2,
            bond_type: bond_type,
            smallest_ring_size: smallest_ring_size,
        }
    }
}
