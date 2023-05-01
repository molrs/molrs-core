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

impl BondType {
    pub fn from_char(bond_type: char) -> Result<BondType, String> {
        if bond_type == ' ' {
            Ok(BondType::Default)
        } else if bond_type == '-' {
            Ok(BondType::Single)
        } else if bond_type == '/' {
            Ok(BondType::Up)
        } else if bond_type == '\\' {
            Ok(BondType::Down)
        } else if bond_type == ':' {
            Ok(BondType::Aromatic)
        } else if bond_type == '=' {
            Ok(BondType::Double)
        } else if bond_type == '#' {
            Ok(BondType::Triple)
        } else if bond_type == '$' {
            Ok(BondType::Quadruple)
        } else {
            Err(format!("invalid bond char {}", bond_type))
        }
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
