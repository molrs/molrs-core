#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum AtomicSymbol {
    H ,                                                                 He,
    Li, Be,                                         B , C , N , O , F , Ne,
    Na, Mg,                                         Al, Si, P , S , Cl, Ar,
    K , Ca, Sc, Ti, V , Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr, Y , Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I , Xe,
    Cs, Ba,     Hf, Ta, W , Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra,     Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,

            La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
            Ac, Th, Pa, U , Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
}

impl AtomicSymbol {
    pub fn new(symbol: &str) -> Result<AtomicSymbol, String> {
        let symbol = symbol.to_lowercase();
        if symbol == "h" { Ok(AtomicSymbol::H) }
        else if symbol == "he" { Ok(AtomicSymbol::He) }
        else if symbol == "li" { Ok(AtomicSymbol::Li) }
        else if symbol == "be" { Ok(AtomicSymbol::Be) }
        else if symbol == "b" { Ok(AtomicSymbol::B) }
        else if symbol == "c" { Ok(AtomicSymbol::C) }
        else if symbol == "n" { Ok(AtomicSymbol::N) }
        else if symbol == "o" { Ok(AtomicSymbol::O) }
        else if symbol == "f" { Ok(AtomicSymbol::F) }
        else if symbol == "ne" { Ok(AtomicSymbol::Ne) }
        else if symbol == "na" { Ok(AtomicSymbol::Na) }
        else if symbol == "mg" { Ok(AtomicSymbol::Mg) }
        else if symbol == "al" { Ok(AtomicSymbol::Al) }
        else if symbol == "si" { Ok(AtomicSymbol::Si) }
        else if symbol == "p" { Ok(AtomicSymbol::P) }
        else if symbol == "s" { Ok(AtomicSymbol::S) }
        else if symbol == "cl" { Ok(AtomicSymbol::Cl) }
        else if symbol == "ar" { Ok(AtomicSymbol::Ar) }
        else if symbol == "k" { Ok(AtomicSymbol::K) }
        else if symbol == "ca" { Ok(AtomicSymbol::Ca) }
        else if symbol == "sc" { Ok(AtomicSymbol::Sc) }
        else if symbol == "ti" { Ok(AtomicSymbol::Ti) }
        else if symbol == "v" { Ok(AtomicSymbol::V) }
        else if symbol == "cr" { Ok(AtomicSymbol::Cr) }
        else if symbol == "mn" { Ok(AtomicSymbol::Mn) }
        else if symbol == "fe" { Ok(AtomicSymbol::Fe) }
        else if symbol == "co" { Ok(AtomicSymbol::Co) }
        else if symbol == "ni" { Ok(AtomicSymbol::Ni) }
        else if symbol == "cu" { Ok(AtomicSymbol::Cu) }
        else if symbol == "zn" { Ok(AtomicSymbol::Zn) }
        else if symbol == "ga" { Ok(AtomicSymbol::Ga) }
        else if symbol == "ge" { Ok(AtomicSymbol::Ge) }
        else if symbol == "as" { Ok(AtomicSymbol::As) }
        else if symbol == "se" { Ok(AtomicSymbol::Se) }
        else if symbol == "br" { Ok(AtomicSymbol::Br) }
        else if symbol == "kr" { Ok(AtomicSymbol::Kr) }
        else if symbol == "rb" { Ok(AtomicSymbol::Rb) }
        else if symbol == "sr" { Ok(AtomicSymbol::Sr) }
        else if symbol == "y" { Ok(AtomicSymbol::Y) }
        else if symbol == "zr" { Ok(AtomicSymbol::Zr) }
        else if symbol == "nb" { Ok(AtomicSymbol::Nb) }
        else if symbol == "mo" { Ok(AtomicSymbol::Mo) }
        else if symbol == "tc" { Ok(AtomicSymbol::Tc) }
        else if symbol == "ru" { Ok(AtomicSymbol::Ru) }
        else if symbol == "rh" { Ok(AtomicSymbol::Rh) }
        else if symbol == "pd" { Ok(AtomicSymbol::Pd) }
        else if symbol == "ag" { Ok(AtomicSymbol::Ag) }
        else if symbol == "cd" { Ok(AtomicSymbol::Cd) }
        else if symbol == "in" { Ok(AtomicSymbol::In) }
        else if symbol == "sn" { Ok(AtomicSymbol::Sn) }
        else if symbol == "sb" { Ok(AtomicSymbol::Sb) }
        else if symbol == "te" { Ok(AtomicSymbol::Te) }
        else if symbol == "i" { Ok(AtomicSymbol::I) }
        else if symbol == "xe" { Ok(AtomicSymbol::Xe) }
        else if symbol == "cs" { Ok(AtomicSymbol::Cs) }
        else if symbol == "ba" { Ok(AtomicSymbol::Ba) }
        else if symbol == "hf" { Ok(AtomicSymbol::Hf) }
        else if symbol == "ta" { Ok(AtomicSymbol::Ta) }
        else if symbol == "w" { Ok(AtomicSymbol::W) }
        else if symbol == "re" { Ok(AtomicSymbol::Re) }
        else if symbol == "os" { Ok(AtomicSymbol::Os) }
        else if symbol == "ir" { Ok(AtomicSymbol::Ir) }
        else if symbol == "pt" { Ok(AtomicSymbol::Pt) }
        else if symbol == "au" { Ok(AtomicSymbol::Au) }
        else if symbol == "hg" { Ok(AtomicSymbol::Hg) }
        else if symbol == "tl" { Ok(AtomicSymbol::Tl) }
        else if symbol == "pb" { Ok(AtomicSymbol::Pb) }
        else if symbol == "bi" { Ok(AtomicSymbol::Bi) }
        else if symbol == "po" { Ok(AtomicSymbol::Po) }
        else if symbol == "at" { Ok(AtomicSymbol::At) }
        else if symbol == "rn" { Ok(AtomicSymbol::Rn) }
        else if symbol == "fr" { Ok(AtomicSymbol::Fr) }
        else if symbol == "ra" { Ok(AtomicSymbol::Ra) }
        else if symbol == "rf" { Ok(AtomicSymbol::Rf) }
        else if symbol == "db" { Ok(AtomicSymbol::Db) }
        else if symbol == "sg" { Ok(AtomicSymbol::Sg) }
        else if symbol == "bh" { Ok(AtomicSymbol::Bh) }
        else if symbol == "hs" { Ok(AtomicSymbol::Hs) }
        else if symbol == "mt" { Ok(AtomicSymbol::Mt) }
        else if symbol == "ds" { Ok(AtomicSymbol::Ds) }
        else if symbol == "rg" { Ok(AtomicSymbol::Rg) }
        else if symbol == "cn" { Ok(AtomicSymbol::Cn) }
        else if symbol == "nh" { Ok(AtomicSymbol::Nh) }
        else if symbol == "fl" { Ok(AtomicSymbol::Fl) }
        else if symbol == "mc" { Ok(AtomicSymbol::Mc) }
        else if symbol == "lv" { Ok(AtomicSymbol::Lv) }
        else if symbol == "ts" { Ok(AtomicSymbol::Ts) }
        else if symbol == "og" { Ok(AtomicSymbol::Og) }
        else if symbol == "la" { Ok(AtomicSymbol::La) }
        else if symbol == "ce" { Ok(AtomicSymbol::Ce) }
        else if symbol == "pr" { Ok(AtomicSymbol::Pr) }
        else if symbol == "nd" { Ok(AtomicSymbol::Nd) }
        else if symbol == "pm" { Ok(AtomicSymbol::Pm) }
        else if symbol == "sm" { Ok(AtomicSymbol::Sm) }
        else if symbol == "eu" { Ok(AtomicSymbol::Eu) }
        else if symbol == "gd" { Ok(AtomicSymbol::Gd) }
        else if symbol == "tb" { Ok(AtomicSymbol::Tb) }
        else if symbol == "dy" { Ok(AtomicSymbol::Dy) }
        else if symbol == "ho" { Ok(AtomicSymbol::Ho) }
        else if symbol == "er" { Ok(AtomicSymbol::Er) }
        else if symbol == "tm" { Ok(AtomicSymbol::Tm) }
        else if symbol == "yb" { Ok(AtomicSymbol::Yb) }
        else if symbol == "lu" { Ok(AtomicSymbol::Lu) }
        else if symbol == "ac" { Ok(AtomicSymbol::Ac) }
        else if symbol == "th" { Ok(AtomicSymbol::Th) }
        else if symbol == "pa" { Ok(AtomicSymbol::Pa) }
        else if symbol == "u" { Ok(AtomicSymbol::U) }
        else if symbol == "np" { Ok(AtomicSymbol::Np) }
        else if symbol == "pu" { Ok(AtomicSymbol::Pu) }
        else if symbol == "am" { Ok(AtomicSymbol::Am) }
        else if symbol == "cm" { Ok(AtomicSymbol::Cm) }
        else if symbol == "bk" { Ok(AtomicSymbol::Bk) }
        else if symbol == "cf" { Ok(AtomicSymbol::Cf) }
        else if symbol == "es" { Ok(AtomicSymbol::Es) }
        else if symbol == "fm" { Ok(AtomicSymbol::Fm) }
        else if symbol == "md" { Ok(AtomicSymbol::Md) }
        else if symbol == "no" { Ok(AtomicSymbol::No) }
        else if symbol == "lr" { Ok(AtomicSymbol::Lr) }
        else { Err(symbol.to_string()) }
    }

    pub fn to_string(&self) -> &str {
        match self {
            AtomicSymbol::H => "H",
            AtomicSymbol::He => "He",
            AtomicSymbol::Li => "Li",
            AtomicSymbol::Be => "Be",
            AtomicSymbol::B => "B",
            AtomicSymbol::C => "C",
            AtomicSymbol::N => "N",
            AtomicSymbol::O => "O",
            AtomicSymbol::F => "F",
            AtomicSymbol::Ne => "Ne",
            AtomicSymbol::Na => "Na",
            AtomicSymbol::Mg => "Mg",
            AtomicSymbol::Al => "Al",
            AtomicSymbol::Si => "Si",
            AtomicSymbol::P => "P",
            AtomicSymbol::S => "S",
            AtomicSymbol::Cl => "Cl",
            AtomicSymbol::Ar => "Ar",
            AtomicSymbol::K => "K",
            AtomicSymbol::Ca => "Ca",
            AtomicSymbol::Sc => "Sc",
            AtomicSymbol::Ti => "Ti",
            AtomicSymbol::V => "V",
            AtomicSymbol::Cr => "Cr",
            AtomicSymbol::Mn => "Mn",
            AtomicSymbol::Fe => "Fe",
            AtomicSymbol::Co => "Co",
            AtomicSymbol::Ni => "Ni",
            AtomicSymbol::Cu => "Cu",
            AtomicSymbol::Zn => "Zn",
            AtomicSymbol::Ga => "Ga",
            AtomicSymbol::Ge => "Ge",
            AtomicSymbol::As => "As",
            AtomicSymbol::Se => "Se",
            AtomicSymbol::Br => "Br",
            AtomicSymbol::Kr => "Kr",
            AtomicSymbol::Rb => "Rb",
            AtomicSymbol::Sr => "Sr",
            AtomicSymbol::Y => "Y",
            AtomicSymbol::Zr => "Zr",
            AtomicSymbol::Nb => "Nb",
            AtomicSymbol::Mo => "Mo",
            AtomicSymbol::Tc => "Tc",
            AtomicSymbol::Ru => "Ru",
            AtomicSymbol::Rh => "Rh",
            AtomicSymbol::Pd => "Pd",
            AtomicSymbol::Ag => "Ag",
            AtomicSymbol::Cd => "Cd",
            AtomicSymbol::In => "In",
            AtomicSymbol::Sn => "Sn",
            AtomicSymbol::Sb => "Sb",
            AtomicSymbol::Te => "Te",
            AtomicSymbol::I => "I",
            AtomicSymbol::Xe => "Xe",
            AtomicSymbol::Cs => "Cs",
            AtomicSymbol::Ba => "Ba",
            AtomicSymbol::Hf => "Hf",
            AtomicSymbol::Ta => "Ta",
            AtomicSymbol::W => "W",
            AtomicSymbol::Re => "Re",
            AtomicSymbol::Os => "Os",
            AtomicSymbol::Ir => "Ir",
            AtomicSymbol::Pt => "Pt",
            AtomicSymbol::Au => "Au",
            AtomicSymbol::Hg => "Hg",
            AtomicSymbol::Tl => "Tl",
            AtomicSymbol::Pb => "Pb",
            AtomicSymbol::Bi => "Bi",
            AtomicSymbol::Po => "Po",
            AtomicSymbol::At => "At",
            AtomicSymbol::Rn => "Rn",
            AtomicSymbol::Fr => "Fr",
            AtomicSymbol::Ra => "Ra",
            AtomicSymbol::Rf => "Rf",
            AtomicSymbol::Db => "Db",
            AtomicSymbol::Sg => "Sg",
            AtomicSymbol::Bh => "Bh",
            AtomicSymbol::Hs => "Hs",
            AtomicSymbol::Mt => "Mt",
            AtomicSymbol::Ds => "Ds",
            AtomicSymbol::Rg => "Rg",
            AtomicSymbol::Cn => "Cn",
            AtomicSymbol::Nh => "Nh",
            AtomicSymbol::Fl => "Fl",
            AtomicSymbol::Mc => "Mc",
            AtomicSymbol::Lv => "Lv",
            AtomicSymbol::Ts => "Ts",
            AtomicSymbol::Og => "Og",
            AtomicSymbol::La => "La",
            AtomicSymbol::Ce => "Ce",
            AtomicSymbol::Pr => "Pr",
            AtomicSymbol::Nd => "Nd",
            AtomicSymbol::Pm => "Pm",
            AtomicSymbol::Sm => "Sm",
            AtomicSymbol::Eu => "Eu",
            AtomicSymbol::Gd => "Gd",
            AtomicSymbol::Tb => "Tb",
            AtomicSymbol::Dy => "Dy",
            AtomicSymbol::Ho => "Ho",
            AtomicSymbol::Er => "Er",
            AtomicSymbol::Tm => "Tm",
            AtomicSymbol::Yb => "Yb",
            AtomicSymbol::Lu => "Lu",
            AtomicSymbol::Ac => "Ac",
            AtomicSymbol::Th => "Th",
            AtomicSymbol::Pa => "Pa",
            AtomicSymbol::U => "U",
            AtomicSymbol::Np => "Np",
            AtomicSymbol::Pu => "Pu",
            AtomicSymbol::Am => "Am",
            AtomicSymbol::Cm => "Cm",
            AtomicSymbol::Bk => "Bk",
            AtomicSymbol::Cf => "Cf",
            AtomicSymbol::Es => "Es",
            AtomicSymbol::Fm => "Fm",
            AtomicSymbol::Md => "Md",
            AtomicSymbol::No => "No",
            AtomicSymbol::Lr => "Lr",
        }
    }

    pub fn atomic_number(&self) -> usize {
        match self {
            AtomicSymbol::H => 1,
            AtomicSymbol::He => 2,
            AtomicSymbol::Li => 3,
            AtomicSymbol::Be => 4,
            AtomicSymbol::B => 5,
            AtomicSymbol::C => 6,
            AtomicSymbol::N => 7,
            AtomicSymbol::O => 8,
            AtomicSymbol::F => 9,
            AtomicSymbol::Ne => 10,
            AtomicSymbol::Na => 11,
            AtomicSymbol::Mg => 12,
            AtomicSymbol::Al => 13,
            AtomicSymbol::Si => 14,
            AtomicSymbol::P => 15,
            AtomicSymbol::S => 16,
            AtomicSymbol::Cl => 17,
            AtomicSymbol::Ar => 18,
            AtomicSymbol::K => 19,
            AtomicSymbol::Ca => 20,
            AtomicSymbol::Sc => 21,
            AtomicSymbol::Ti => 22,
            AtomicSymbol::V => 23,
            AtomicSymbol::Cr => 24,
            AtomicSymbol::Mn => 25,
            AtomicSymbol::Fe => 26,
            AtomicSymbol::Co => 27,
            AtomicSymbol::Ni => 28,
            AtomicSymbol::Cu => 29,
            AtomicSymbol::Zn => 30,
            AtomicSymbol::Ga => 31,
            AtomicSymbol::Ge => 32,
            AtomicSymbol::As => 33,
            AtomicSymbol::Se => 34,
            AtomicSymbol::Br => 35,
            AtomicSymbol::Kr => 36,
            AtomicSymbol::Rb => 37,
            AtomicSymbol::Sr => 38,
            AtomicSymbol::Y => 39,
            AtomicSymbol::Zr => 40,
            AtomicSymbol::Nb => 41,
            AtomicSymbol::Mo => 42,
            AtomicSymbol::Tc => 43,
            AtomicSymbol::Ru => 44,
            AtomicSymbol::Rh => 45,
            AtomicSymbol::Pd => 46,
            AtomicSymbol::Ag => 47,
            AtomicSymbol::Cd => 48,
            AtomicSymbol::In => 49,
            AtomicSymbol::Sn => 50,
            AtomicSymbol::Sb => 51,
            AtomicSymbol::Te => 52,
            AtomicSymbol::I => 53,
            AtomicSymbol::Xe => 54,
            AtomicSymbol::Cs => 55,
            AtomicSymbol::Ba => 56,
            AtomicSymbol::Hf => 72,
            AtomicSymbol::Ta => 73,
            AtomicSymbol::W => 74,
            AtomicSymbol::Re => 75,
            AtomicSymbol::Os => 76,
            AtomicSymbol::Ir => 77,
            AtomicSymbol::Pt => 78,
            AtomicSymbol::Au => 79,
            AtomicSymbol::Hg => 80,
            AtomicSymbol::Tl => 81,
            AtomicSymbol::Pb => 82,
            AtomicSymbol::Bi => 83,
            AtomicSymbol::Po => 84,
            AtomicSymbol::At => 85,
            AtomicSymbol::Rn => 86,
            AtomicSymbol::Fr => 87,
            AtomicSymbol::Ra => 88,
            AtomicSymbol::Rf => 104,
            AtomicSymbol::Db => 105,
            AtomicSymbol::Sg => 106,
            AtomicSymbol::Bh => 107,
            AtomicSymbol::Hs => 108,
            AtomicSymbol::Mt => 109,
            AtomicSymbol::Ds => 110,
            AtomicSymbol::Rg => 111,
            AtomicSymbol::Cn => 112,
            AtomicSymbol::Nh => 113,
            AtomicSymbol::Fl => 114,
            AtomicSymbol::Mc => 115,
            AtomicSymbol::Lv => 116,
            AtomicSymbol::Ts => 117,
            AtomicSymbol::Og => 118,
            AtomicSymbol::La => 57,
            AtomicSymbol::Ce => 58,
            AtomicSymbol::Pr => 59,
            AtomicSymbol::Nd => 60,
            AtomicSymbol::Pm => 61,
            AtomicSymbol::Sm => 62,
            AtomicSymbol::Eu => 63,
            AtomicSymbol::Gd => 64,
            AtomicSymbol::Tb => 65,
            AtomicSymbol::Dy => 66,
            AtomicSymbol::Ho => 67,
            AtomicSymbol::Er => 68,
            AtomicSymbol::Tm => 69,
            AtomicSymbol::Yb => 70,
            AtomicSymbol::Lu => 71,
            AtomicSymbol::Ac => 89,
            AtomicSymbol::Th => 90,
            AtomicSymbol::Pa => 91,
            AtomicSymbol::U => 92,
            AtomicSymbol::Np => 93,
            AtomicSymbol::Pu => 94,
            AtomicSymbol::Am => 95,
            AtomicSymbol::Cm => 96,
            AtomicSymbol::Bk => 97,
            AtomicSymbol::Cf => 98,
            AtomicSymbol::Es => 99,
            AtomicSymbol::Fm => 100,
            AtomicSymbol::Md => 101,
            AtomicSymbol::No => 102,
            AtomicSymbol::Lr => 103,
        }
    }

    pub fn default_isotope(&self) -> usize {
        match self {
            AtomicSymbol::H => 1,
            AtomicSymbol::He => 4,
            AtomicSymbol::Li => 7,
            AtomicSymbol::Be => 9,
            AtomicSymbol::B => 11,
            AtomicSymbol::C => 12,
            AtomicSymbol::N => 14,
            AtomicSymbol::O => 16,
            AtomicSymbol::F => 19,
            AtomicSymbol::Ne => 20,
            AtomicSymbol::Na => 23,
            AtomicSymbol::Mg => 24,
            AtomicSymbol::Al => 27,
            AtomicSymbol::Si => 28,
            AtomicSymbol::P => 31,
            AtomicSymbol::S => 32,
            AtomicSymbol::Cl => 35,
            AtomicSymbol::Ar => 40,
            AtomicSymbol::K => 39,
            AtomicSymbol::Ca => 40,
            AtomicSymbol::Sc => 45,
            AtomicSymbol::Ti => 48,
            AtomicSymbol::V => 51,
            AtomicSymbol::Cr => 52,
            AtomicSymbol::Mn => 55,
            AtomicSymbol::Fe => 56,
            AtomicSymbol::Co => 59,
            AtomicSymbol::Ni => 59,
            AtomicSymbol::Cu => 64,
            AtomicSymbol::Zn => 65,
            AtomicSymbol::Ga => 70,
            AtomicSymbol::Ge => 73,
            AtomicSymbol::As => 75,
            AtomicSymbol::Se => 79,
            AtomicSymbol::Br => 80,
            AtomicSymbol::Kr => 84,
            AtomicSymbol::Rb => 85,
            AtomicSymbol::Sr => 88,
            AtomicSymbol::Y => 89,
            AtomicSymbol::Zr => 91,
            AtomicSymbol::Nb => 93,
            AtomicSymbol::Mo => 96,
            AtomicSymbol::Tc => 98,
            AtomicSymbol::Ru => 101,
            AtomicSymbol::Rh => 103,
            AtomicSymbol::Pd => 106,
            AtomicSymbol::Ag => 108,
            AtomicSymbol::Cd => 112,
            AtomicSymbol::In => 115,
            AtomicSymbol::Sn => 119,
            AtomicSymbol::Sb => 122,
            AtomicSymbol::Te => 128,
            AtomicSymbol::I => 127,
            AtomicSymbol::Xe => 131,
            AtomicSymbol::Cs => 133,
            AtomicSymbol::Ba => 137,
            AtomicSymbol::Hf => 178,
            AtomicSymbol::Ta => 181,
            AtomicSymbol::W => 184,
            AtomicSymbol::Re => 186,
            AtomicSymbol::Os => 190,
            AtomicSymbol::Ir => 192,
            AtomicSymbol::Pt => 195,
            AtomicSymbol::Au => 197,
            AtomicSymbol::Hg => 201,
            AtomicSymbol::Tl => 204,
            AtomicSymbol::Pb => 207,
            AtomicSymbol::Bi => 209,
            AtomicSymbol::Po => 209,
            AtomicSymbol::At => 210,
            AtomicSymbol::Rn => 222,
            AtomicSymbol::Fr => 223,
            AtomicSymbol::Ra => 226,
            AtomicSymbol::Rf => 267,
            AtomicSymbol::Db => 268,
            AtomicSymbol::Sg => 269,
            AtomicSymbol::Bh => 270,
            AtomicSymbol::Hs => 269,
            AtomicSymbol::Mt => 278,
            AtomicSymbol::Ds => 281,
            AtomicSymbol::Rg => 281,
            AtomicSymbol::Cn => 285,
            AtomicSymbol::Nh => 284,
            AtomicSymbol::Fl => 289,
            AtomicSymbol::Mc => 288,
            AtomicSymbol::Lv => 293,
            AtomicSymbol::Ts => 292,
            AtomicSymbol::Og => 294,
            AtomicSymbol::La => 139,
            AtomicSymbol::Ce => 140,
            AtomicSymbol::Pr => 141,
            AtomicSymbol::Nd => 144,
            AtomicSymbol::Pm => 145,
            AtomicSymbol::Sm => 150,
            AtomicSymbol::Eu => 152,
            AtomicSymbol::Gd => 157,
            AtomicSymbol::Tb => 159,
            AtomicSymbol::Dy => 162,
            AtomicSymbol::Ho => 165,
            AtomicSymbol::Er => 167,
            AtomicSymbol::Tm => 169,
            AtomicSymbol::Yb => 173,
            AtomicSymbol::Lu => 175,
            AtomicSymbol::Ac => 227,
            AtomicSymbol::Th => 232,
            AtomicSymbol::Pa => 231,
            AtomicSymbol::U => 238,
            AtomicSymbol::Np => 237,
            AtomicSymbol::Pu => 244,
            AtomicSymbol::Am => 243,
            AtomicSymbol::Cm => 247,
            AtomicSymbol::Bk => 247,
            AtomicSymbol::Cf => 251,
            AtomicSymbol::Es => 252,
            AtomicSymbol::Fm => 257,
            AtomicSymbol::Md => 258,
            AtomicSymbol::No => 259,
            AtomicSymbol::Lr => 262,
        }
    }

    pub fn atomic_mass(&self) {}

    pub fn num_imp_h(&self) -> usize {
        match self {
            AtomicSymbol::B => 3,
            AtomicSymbol::C => 4,
            AtomicSymbol::N => 3,
            AtomicSymbol::O => 2,
            AtomicSymbol::F => 1,
            AtomicSymbol::P => 3,
            AtomicSymbol::S => 2,
            AtomicSymbol::Cl => 1,
            AtomicSymbol::Br => 1,
            AtomicSymbol::I => 1,
            _ => 0,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn atomic_symbol_from_str() {
        assert_eq!(AtomicSymbol::new("H").unwrap(), AtomicSymbol::H);
        assert_eq!(AtomicSymbol::new("h").unwrap(), AtomicSymbol::H);
        assert_eq!(AtomicSymbol::new("He").unwrap(), AtomicSymbol::He);
        assert_eq!(AtomicSymbol::new("he").unwrap(), AtomicSymbol::He);
    }
    #[test]
    fn atomic_symbol_to_str() {
        assert_eq!(AtomicSymbol::new("H").unwrap().to_string(), "H");
        assert_eq!(AtomicSymbol::new("he").unwrap().to_string(), "He");
    }
    #[test]
    fn atomic_symbol_to_atomic_number() {
        assert_eq!(AtomicSymbol::new("H").unwrap().atomic_number(), 1);
        assert_eq!(AtomicSymbol::new("he").unwrap().atomic_number(), 2);
    }
}
