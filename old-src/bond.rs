use std::cmp::Ordering;

#[derive(Debug)]
pub struct BondParseError {
    pub details: String,
}

/// Enum that represents the bond type.
#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum BondType {
    #[default]
    Default,
    Single,
    Up,
    Down,
    Delocalized,
    Double,
    Triple,
    Quadruple,
    Dative,
}

impl TryFrom<char> for BondType {
    type Error = BondParseError;

    /// Matches chars to BondTypes following the OpenSMILES syntax.
    fn try_from(bond_type: char) -> Result<Self, Self::Error> {
        match bond_type {
            ' ' => Ok(BondType::Default),
            '-' => Ok(BondType::Single),
            '/' => Ok(BondType::Up),
            '\\' => Ok(BondType::Down),
            ':' => Ok(BondType::Delocalized),
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
    /// Returns the bond order.
    pub fn to_float(&self) -> f64 {
        match self {
            BondType::Default => 1.0,
            BondType::Single => 1.0,
            BondType::Up => 1.0,
            BondType::Down => 1.0,
            BondType::Delocalized => 1.5,
            BondType::Double => 2.0,
            BondType::Triple => 3.0,
            BondType::Quadruple => 4.0,
            BondType::Dative => unimplemented!(),
        }
    }

    /// Returns the bond type as a SMILES char.
    pub fn to_char(&self) -> char {
        match self {
            BondType::Default => ' ',
            BondType::Single => '-',
            BondType::Up => '/',
            BondType::Down => '\\',
            BondType::Delocalized => ':',
            BondType::Double => '=',
            BondType::Triple => '#',
            BondType::Quadruple => '$',
            BondType::Dative => unimplemented!(),
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Bond {
    pub atom_i: usize,
    pub atom_j: usize,
    pub bond_type: BondType,
}

impl Ord for Bond {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.atom_i.cmp(&other.atom_i) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => match self.atom_j.cmp(&other.atom_j) {
                Ordering::Greater => Ordering::Greater,
                Ordering::Less => Ordering::Less,
                Ordering::Equal => Ordering::Equal,
            },
        }
    }
}

impl PartialOrd for Bond {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Bond {
    pub fn new(atom_i: usize, atom_j: usize, bond_type: char) -> Result<Bond, BondParseError> {
        let bond_type = BondType::try_from(bond_type)?;
        match atom_i.cmp(&atom_j) {
            Ordering::Greater => Ok(Bond {
                atom_j,
                atom_i,
                bond_type,
            }),
            Ordering::Less => Ok(Bond {
                atom_i,
                atom_j,
                bond_type,
            }),
            Ordering::Equal => Err(BondParseError {
                details: "atom cannot bond to itself".to_owned(),
            }),
        }
    }
}
