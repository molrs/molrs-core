use pertable::Element;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PointChirality {
    #[default]
    Undefined,
    Achiral,
    Clockwise,
    CounterClockwise,
}

impl ToString for PointChirality {
    fn to_string(&self) -> String {
        match self {
            PointChirality::Undefined => "",
            PointChirality::Achiral => "",
            PointChirality::CounterClockwise => "@",
            PointChirality::Clockwise => "@@",
        }
        .to_owned()
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

impl ToString for Atom {
    fn to_string(&self) -> String {
        let mut s = String::new();
        if self.isotope.is_none()
            && self.charge == 0
            && self.n_radical_electrons == Some(0)
            && (self.point_chirality == PointChirality::Undefined
                || self.point_chirality == PointChirality::Achiral)
        {
            let atomic_symbol = match self.delocalized {
                true => self.element.to_string().to_lowercase(),
                false => self.element.to_string(),
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

        s
    }
}
