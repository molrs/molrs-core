use pertable::Element;
use std::fmt;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PointChirality {
    #[default]
    Undefined,
    Achiral,
    Clockwise,
    CounterClockwise,
}

impl fmt::Display for PointChirality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            PointChirality::Undefined => "",
            PointChirality::Achiral => "",
            PointChirality::CounterClockwise => "@",
            PointChirality::Clockwise => "@@",
        };
        write!(f, "{s}",)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Atom {
    pub element: Element,
    pub isotope: Option<u16>,
    pub charge: i8,
    pub delocalized: bool,
    pub n_implicit_hydrogens: Option<u8>,
    pub n_radical_electrons: Option<u8>,
    pub point_chirality: PointChirality,
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        if self.isotope.is_none()
            && self.charge == 0
            && (self.n_radical_electrons == Some(0) || self.n_radical_electrons.is_none())
            && (self.point_chirality == PointChirality::Undefined
                || self.point_chirality == PointChirality::Achiral)
        {
            let atomic_symbol = match self.delocalized {
                true => self.element.atomic_symbol().to_lowercase(),
                false => self.element.atomic_symbol(),
            };

            s = match self.element {
                Element::Any
                | Element::B
                | Element::C
                | Element::N
                | Element::O
                | Element::S
                | Element::P
                | Element::F
                | Element::Cl
                | Element::Br
                | Element::I => atomic_symbol,
                _ => format!("[{atomic_symbol}]"),
            };
        } else {
            s.push('[');
            if self.isotope.is_some() {
                s += &self.isotope.unwrap().to_string();
            };
            s += &self.element.to_string();
            if self.point_chirality != PointChirality::Undefined
                || self.point_chirality != PointChirality::Achiral
            {
                s += &self.point_chirality.to_string();
            }
            s += &match self.n_implicit_hydrogens {
                None => "".to_owned(),
                Some(n_implicit_hydrogens) => match n_implicit_hydrogens {
                    0 => "".to_owned(),
                    1 => "H".to_owned(),
                    _ => format!("H{n_implicit_hydrogens}"),
                },
            };
            if self.charge == 1 {
                s += "+";
            } else if self.charge == -1 {
                s += "-";
            } else if self.charge > 0 {
                s += &format!("+{}", self.charge);
            } else if self.charge < 0 {
                s += &format!("-{}", self.charge);
            }
            s.push(']');
        }

        write! {f, "{s}"}
    }
}
