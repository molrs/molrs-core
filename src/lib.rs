#![feature(test)]

pub mod atom;
pub mod extended_atom;
pub mod from;
pub mod molecule;
pub mod reaction;
pub mod read;
pub mod standardize;
pub mod substructure;
pub mod to;
pub mod write;

pub fn visualize() {}

#[cfg(test)]
mod tests {
    extern crate test;
    use super::*;

    #[test]
    fn test() {
        println!("meep");
    }
}
