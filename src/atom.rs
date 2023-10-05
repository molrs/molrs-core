use pertable::Element;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PointChirality {
    #[default]
    Undefined,
    Achiral,
    Clockwise,
    CounterClockwise,
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
