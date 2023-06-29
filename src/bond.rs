#[derive(Debug)]
pub struct BondParseError {
    pub details: String,
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum BondType {
    #[default]
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

impl TryFrom<char> for BondType {
    type Error = BondParseError;

    fn try_from(bond_type: char) -> Result<Self, Self::Error> {
        match bond_type {
            ' ' => Ok(BondType::Default),
            '-' => Ok(BondType::Single),
            '/' => Ok(BondType::Up),
            '\\' => Ok(BondType::Down),
            ':' => Ok(BondType::Aromatic),
            '=' => Ok(BondType::Double),
            '#' => Ok(BondType::Triple),
            '$' => Ok(BondType::Quadruple),
            _ => Err(BondParseError {
                details: format!("invalid bond char {}", bond_type),
            }),
        }
    }
}

impl BondType {
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
    atom_i: usize,
    atom_j: usize,
    bond_type: BondType,
}

impl Bond {
    pub fn new(atom_i: usize, atom_j: usize, bond_type: char) -> Result<Bond, BondParseError> {
        let bond_type = BondType::try_from(bond_type)?;
        Ok(Bond {
            atom_i,
            atom_j,
            bond_type,
        })
    }
}
