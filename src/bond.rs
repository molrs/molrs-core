use std::cmp::Ordering;

#[derive(Debug)]
pub enum BondError {
    InvalidBondChar(char),
    AtomBondedToSelf,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
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
    type Error = BondError;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            ' ' => Ok(BondType::Default),
            '-' => Ok(BondType::Single),
            '/' => Ok(BondType::Up),
            '\\' => Ok(BondType::Down),
            ':' => Ok(BondType::Delocalized),
            '=' => Ok(BondType::Double),
            '#' => Ok(BondType::Triple),
            '$' => Ok(BondType::Quadruple),
            _ => Err(BondError::InvalidBondChar(value)),
        }
    }
}

impl From<BondType> for char {
    fn from(value: BondType) -> Self {
        match value {
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

impl From<&BondType> for char {
    fn from(value: &BondType) -> Self {
        match value {
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

impl From<BondType> for f64 {
    fn from(value: BondType) -> Self {
        match value {
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
}

impl From<&BondType> for f64 {
    fn from(value: &BondType) -> Self {
        match value {
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
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub bond_type: BondType,
}

impl Bond {
    pub fn new(i: usize, j: usize, bond_type: char) -> Result<Bond, BondError> {
        let bond_type = BondType::try_from(bond_type)?;
        match i.cmp(&j) {
            Ordering::Greater => Ok(Bond { j, i, bond_type }),
            Ordering::Less => Ok(Bond { i, j, bond_type }),
            Ordering::Equal => Err(BondError::AtomBondedToSelf),
        }
    }
}
