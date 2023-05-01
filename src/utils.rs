use std::{collections::HashSet};

pub fn get_duplicate_element(vector: &Vec<usize>) -> Option<usize> {
    let mut set = HashSet::new();
    for value in vector {
        if set.contains(value) {
            return Some(*value);
        } else {
            set.insert(value);
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_duplicate_elements() {
        assert_eq!(get_duplicate_element(&vec![0, 1, 2, 3, 4, 0]).unwrap(), 0);
        assert_eq!(get_duplicate_element(&vec![0, 1, 2, 3, 4]), Option::None);
    }
}
