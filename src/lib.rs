#![feature(test)]
#![feature(iterator_try_collect)]

pub mod atom;
pub mod bond;
pub mod element;
pub mod molecule;
pub mod smiles;
pub mod utils;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let e = element::Element::Ac;
        let e = String::from(&e);
        dbg!(&e);
    }
}
