#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum BondType {
    Default,
    Single,
    Up,
    Down,
    Aromatic,
    Double,
    Triple,
    Quadruple,
    Dative,
}

impl BondType {
    pub fn new(bond_type: char) -> Result<BondType, String> {
        if bond_type == ' ' { Ok(BondType::Default) }
        else if bond_type == '-' { Ok(BondType::Single) }
        else if bond_type == '/' { Ok(BondType::Up) }
        else if bond_type == '\\' { Ok(BondType::Down) }
        else if bond_type == ':' { Ok(BondType::Aromatic) }
        else if bond_type == '=' { Ok(BondType::Double) }
        else if bond_type == '#' { Ok(BondType::Triple) }
        else if bond_type == '$' { Ok(BondType::Quadruple) }
        else if bond_type == '~' { Ok(BondType::Dative) }
        else { Err(format!("invalid bond_type, {}", bond_type)) }
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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Bond {
    pub atom_idx_1: usize,
    pub atom_idx_2: usize,
    pub bond_type: BondType,
}

impl Bond {
    pub fn new(atom_idx_1: usize, atom_idx_2: usize, bond_type: BondType) -> Bond {
        Bond {
            atom_idx_1: atom_idx_1,
            atom_idx_2: atom_idx_2,
            bond_type: bond_type,
        }
    }

    pub fn from_char(atom_idx_1: usize, atom_idx_2: usize, bond_type: char) -> Result<Bond, String> {
        let bond_type = match BondType::new(bond_type) {
            Ok(link_type) => link_type,
            Err(string) => return Err(string),
        };
        Ok(Bond::new(atom_idx_1, atom_idx_2, bond_type))
    }
}
