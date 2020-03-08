use std::hint::unreachable_unchecked;
use std::collections::{BTreeMap, BTreeSet};
use std::num::NonZeroU8;

use ptable::Element;

use super::*;
use crate::parse_number;

#[derive(Debug)]
pub struct SMILESMisc {
    pub automatic_hydrogens_targets: Vec<usize>,
    pub labels: BTreeMap<u32, usize>,
    pub aromatic_ions: BTreeSet<usize>,
}

#[derive(Debug, Clone)]
pub struct AromaticDetectionData {
    pub init: bool,
    pub last_one_was_double: bool,
}

// aromatic_detection_data = Small structure composed of information to the formation of an aromatic ring,
// bonds_in_the_probably_aromatic_ring = The atoms that possibly make up the aromatic ring that is being checked, this is restored by the function called
#[allow(unused_assignments)]  // Rust compiler is buggy
fn parse_smiles_group<'a>(mut string: &'a [u8], sf: &mut Molecule, misc: &mut SMILESMisc, 
    mut adj_atom: Option<usize>, aromatic_detection_data: &mut AromaticDetectionData,
    bonds_in_the_probably_aromatic_ring: &mut Vec<usize>) -> Result<&'a [u8], ()> {
    const SMILES_SPECIAL_UNICHR_ELEMENTS: [(u8, Element); 6] = [(b'F', Element::Fluorine), (b'I', Element::Iodine), 
        (b'N', Element::Nitrogen), (b'O', Element::Oxygen), (b'P', Element::Phosphorus), (b'S', Element::Sulfur)];
    let mut waiting_bond = None; // The bond is parse first and store here for the next ion, which takes it
    let bond_count_in_the_probably_aromatic_ring = bonds_in_the_probably_aromatic_ring.len();
    macro_rules! insert_ion {
        ($atom:expr, $is_aromatic:expr, $charge:expr, $isotopic_spec:expr) => {{
            let mut labels = Vec::new();
            while !string.is_empty() {
                let mut bond = None;
                let mut pos = 0;
                if string[0] == b':' || string[0] == b'-' || string[0] == b'=' || string[0] == b'#' {
                    bond = Some(match string[0] {
                        b':' => StructuralBond::Aromatic,
                        b'-' => StructuralBond::Single,
                        b'=' => StructuralBond::Double,
                        b'#' => StructuralBond::Triple,
                        _ => unsafe { unreachable_unchecked() }
                    });
                    pos += 1;
                    if string.len() == 1 {
                        return Err(());
                    }
                }
                if string[pos] >= b'0' && string[pos] <= b'9' {
                    // One-digit label
                    labels.push(((string[pos] - b'0') as u32, bond));
                    string = &string[pos+1..];
                }
                else if string[pos] == b'%' {
                    // Multi-digit label
                    string = &string[pos+1..];
                    if string.is_empty() || string[0] < b'0' || string[0] > b'9' {
                        return Err(());
                    }
                    let (n, s) = parse_number(string);
                    string = s;
                    labels.push((n as u32, bond));
                }
                else {
                    break;
                }
            }
            let mut current_bounds = Vec::new();
            if let Some(adj_atom) = adj_atom {
                let k; // k = Bound kind if this bond
                if let Some(bond) = waiting_bond {
                    k = bond;
                    waiting_bond = None;
                }
                else {
                    if $is_aromatic && misc.aromatic_ions.contains(&adj_atom) {
                        k = StructuralBond::Aromatic;
                    }
                    else {
                        k = StructuralBond::Single;
                    }
                }
                current_bounds.push(sf.bonds.len());
                sf.atoms[adj_atom].bonds.push(sf.bonds.len());
                if aromatic_detection_data.init {
                    // Default believe there is no aromaticity
                    aromatic_detection_data.init = false;
                    if k == StructuralBond::Single || k == StructuralBond::Double {
                        if ((k == StructuralBond::Single) == aromatic_detection_data.last_one_was_double) {
                            // =X-C or -X=C bond
                            // There is aromaticity, so far
                            aromatic_detection_data.last_one_was_double = !aromatic_detection_data.last_one_was_double;
                            aromatic_detection_data.init = true;
                            bonds_in_the_probably_aromatic_ring.push(sf.bonds.len());
                        }
                    }
                    if !aromatic_detection_data.init {
                        bonds_in_the_probably_aromatic_ring.truncate(bond_count_in_the_probably_aromatic_ring);
                    }
                }
                else {
                    if k == StructuralBond::Single || k == StructuralBond::Double {
                        aromatic_detection_data.init = true;
                        aromatic_detection_data.last_one_was_double = k == StructuralBond::Double;
                        bonds_in_the_probably_aromatic_ring.push(sf.bonds.len());
                    }
                }
                sf.bonds.push(Bond { a: adj_atom, b: sf.atoms.len(), k });
                waiting_bond = None;
            }
            for label in labels {
                if let Some(b) = misc.labels.get(&label.0) {
                    let b = *b;
                    let mut k = match label.1 {
                        Some(k) => k,
                        None => {
                            if $is_aromatic && misc.aromatic_ions.contains(&b) {
                                StructuralBond::Aromatic
                            }
                            else {
                                StructuralBond::Single
                            }
                        }
                    };
                    misc.labels.remove(&label.0);
                    if aromatic_detection_data.init && bonds_in_the_probably_aromatic_ring.len() >= 4*1+2-1 {
                        // Seek b atom in the possible ring
                        let mut other_end = None;
                        for (p, id) in bonds_in_the_probably_aromatic_ring.iter().enumerate() {
                            if *id == b {
                                other_end = Some(p);
                                break;
                            }
                        }
                        if let Some(other_end) = other_end {
                            let length = bonds_in_the_probably_aromatic_ring.len() - other_end + 1;
                            if (length - 2) % 4 == 0 { // Test if it complies the Hueckel's rule
                                for i in bonds_in_the_probably_aromatic_ring.iter() {
                                    let i = *i;
                                    sf.bonds[i].k = StructuralBond::Aromatic;
                                }
                                k = StructuralBond::Aromatic;
                                bonds_in_the_probably_aromatic_ring.clear();
                                aromatic_detection_data.init = false;
                            }
                        }
                    }
                    current_bounds.push(sf.bonds.len());
                    sf.atoms[b].bonds.push(sf.bonds.len());
                    sf.bonds.push(Bond { a: sf.atoms.len(), b, k });
                }
                else {
                    misc.labels.insert(label.0, sf.atoms.len());
                }
            }
            if $is_aromatic {
                misc.aromatic_ions.insert(sf.atoms.len());
            }
            adj_atom = Some(sf.atoms.len());
            sf.atoms.push(AtomAndBondI { atom: Isotope::new(Ion::new($atom, $charge), $isotopic_spec), bonds: current_bounds });
        }};
    }
    if let Some(_) = adj_atom {
        if !string.is_empty() {
            match string[0] {
                b'-' => {
                    waiting_bond = Some(StructuralBond::Single);
                    string = &string[1..];
                },
                b'=' => {
                    waiting_bond = Some(StructuralBond::Double);
                    string = &string[1..];
                },
                b'#' => {
                    waiting_bond = Some(StructuralBond::Triple);
                    string = &string[1..];
                },
                _ => {}
            }
        }
    }
    'mainloop: while !string.is_empty() {
        let chr = string[0];
        if chr == b')' {
            break;
        }
        string = &string[1..];
        if chr == b'B' {
            if string.is_empty() {
                misc.automatic_hydrogens_targets.push(sf.atoms.len());
                insert_ion!(Element::Boron, false, 0, None);
                break;
            }
            if string[0] == b'r' {
                string = &string[1..];
                misc.automatic_hydrogens_targets.push(sf.atoms.len());
                insert_ion!(Element::Bromine, false, 0, None);
            }
            else {
                misc.automatic_hydrogens_targets.push(sf.atoms.len());
                insert_ion!(Element::Boron, false, 0, None);
            }
            continue;
        }
        else if chr == b'C' {
            if string.is_empty() {
                misc.automatic_hydrogens_targets.push(sf.atoms.len());
                insert_ion!(Element::Carbon, false, 0, None);
                break;
            }
            if string[0] == b'l' {
                string = &string[1..];
                misc.automatic_hydrogens_targets.push(sf.atoms.len());
                insert_ion!(Element::Chlorine, false, 0, None);
            }
            else {
                misc.automatic_hydrogens_targets.push(sf.atoms.len());
                insert_ion!(Element::Carbon, false, 0, None);
            }
            continue;
        }
        else if chr == b'[' {
            if string.is_empty() {
                return Err(());
            }
            let mut isotopic_spec = None;
            if string[0] >= b'0' && string[0] <= b'9' {
                // Isotopic specification
                let (r, s) = parse_number(string);
                string = s;
                if r >= 256 {
                    return Err(());
                }
                isotopic_spec = NonZeroU8::new(r as u8);
                if string.is_empty() {
                    return Err(());
                }
            }
            while string[0] != b']' {
                if string.len() < 2 { // 1: X, 2: ]  "X]"
                    return Err(());
                }
                let e;
                let mut aromatic = false;
                if string[0] == b'c' || string[0] == b'n' || string[0] == b'o' || string[0] == b's' {
                    e = match string[0] {
                        b'c' => Element::Carbon,
                        b'n' => Element::Nitrogen,
                        b'o' => Element::Phosphorus,
                        b's' => Element::Sulfur,
                        _ => unsafe { unreachable_unchecked() }
                    };
                    aromatic = true;
                    string = &string[1..];
                }
                else if string[1] >= b'a' && string[1] <= b'z' {
                    e = match Element::from_symbol(unsafe { mem::transmute(&string[..2]) }) {
                        Some(e) => e,
                        None => return Err(())
                    };
                    string = &string[2..];
                    if string.is_empty() {
                        return Err(());
                    }
                }
                else {
                    e = match Element::from_symbol(unsafe { mem::transmute(&string[..1]) }) {
                        Some(e) => e,
                        None => return Err(())
                    };
                    string = &string[1..];
                }
                let mut charge = 0;
                if string[0] == b'+' || string[0] == b'-' {
                    let negative = string[0] == b'-';
                    string = &string[1..];
                    if string.is_empty() {
                        return Err(());
                    }
                    if string[0] >= b'0' && string[0] <= b'9' {
                        // Specific with number
                        let (n, s) = parse_number(string);
                        string = s;
                        if string.is_empty() {
                            return Err(());
                        }
                        if n > 127 {
                            return Err(());
                        }
                        charge = n as i8;
                    }
                    else {
                        // Specific with symbol count
                        charge += 1;
                        while !string.is_empty() {
                            if string[0] != b'+' && string[0] != b'-' {
                                break;
                            }
                            if (string[0] == b'-') != negative {
                                // Weird charge specification, yield an error
                                return Err(());
                            }
                            charge += 1;
                        }
                        if string.is_empty() {
                            return Err(());
                        }
                    }
                    if negative {
                        charge = -charge;
                    }
                }
                insert_ion!(e, aromatic, charge, isotopic_spec);
                isotopic_spec = None;
            }
            string = &string[1..];
        }
        else if chr >= b'A' && chr <= b'Z' {
            let mut l = 0;
            let mut r = 5;
            while l <= r {
                let m = l + (r - l) / 2;
                if SMILES_SPECIAL_UNICHR_ELEMENTS[m].0 == chr {
                    misc.automatic_hydrogens_targets.push(sf.atoms.len());
                    insert_ion!(SMILES_SPECIAL_UNICHR_ELEMENTS[m].1, false, 0, None);
                    continue 'mainloop;
                }
                if SMILES_SPECIAL_UNICHR_ELEMENTS[m].0 < chr {
                    l = m + 1;
                }
                if SMILES_SPECIAL_UNICHR_ELEMENTS[m].0 > chr {
                    r = m - 1;
                }
            }
        }
        else if chr == b'c' {
            misc.automatic_hydrogens_targets.push(sf.atoms.len());
            insert_ion!(Element::Carbon, true, 0, None);
        }
        else if chr == b'o' {
            misc.automatic_hydrogens_targets.push(sf.atoms.len());
            insert_ion!(Element::Oxygen, true, 0, None);
        }
        else if chr == b's' {
            misc.automatic_hydrogens_targets.push(sf.atoms.len());
            insert_ion!(Element::Sulfur, true, 0, None);
        }
        else if chr == b'n' {
            misc.automatic_hydrogens_targets.push(sf.atoms.len());
            insert_ion!(Element::Nitrogen, true, 0, None);
        }
        else if chr == b'(' {
            // Branch
            let mut new_data = aromatic_detection_data.clone();
            string = parse_smiles_group(string, sf, misc, adj_atom, &mut new_data, bonds_in_the_probably_aromatic_ring)?;
            if bonds_in_the_probably_aromatic_ring.is_empty() {
                // Aromatic ring close or not started
                aromatic_detection_data.init = false;
            }
            // SMILES ends without close branch
            if string.is_empty() {
                return Err(());
            }
            string = &string[1..]
        }
        else {
            if let None = waiting_bond {
                match chr {
                    b':' => {
                        waiting_bond = Some(StructuralBond::Aromatic);
                        continue;
                    },
                    b'-' => {
                        waiting_bond = Some(StructuralBond::Single);
                        continue;
                    },
                    b'=' => {
                        waiting_bond = Some(StructuralBond::Double);
                        continue;
                    },
                    b'#' => {
                        waiting_bond = Some(StructuralBond::Triple);
                        continue;
                    },
                    _ => {}
                }
            }
            return Err(());
        }
    }
    bonds_in_the_probably_aromatic_ring.truncate(bond_count_in_the_probably_aromatic_ring);
    Ok(string)
}

pub fn parse(string: &str) -> Result<Molecule, ()> {
    let mut sf = Molecule { atoms: Vec::new(), bonds: Vec::new(), coords: None };
    let string = string.as_bytes();
    let mut misc = smiles::SMILESMisc { automatic_hydrogens_targets: Vec::new(), 
        labels: BTreeMap::new(), aromatic_ions: BTreeSet::new() };
    let _ = smiles::parse_smiles_group(string, &mut sf, &mut misc, None, &mut smiles::AromaticDetectionData { init: false, 
        last_one_was_double: false }, &mut Vec::new())?;
    for atom in misc.automatic_hydrogens_targets {
        let mut actual_valence = 0;
        let mut aromatic_flip_flop = true;
        for bond in sf.atoms[atom].bonds.iter() {
            if sf.bonds[*bond].k == StructuralBond::Aromatic {
                actual_valence += 1;
                if aromatic_flip_flop {
                    actual_valence += 1;
                }
                aromatic_flip_flop = !aromatic_flip_flop;
            }
            actual_valence += sf.bonds[*bond].k as isize;
        }
        let lack = match sf.atoms[atom].get_element() {
            Element::Boron => 3 - actual_valence,
            Element::Carbon => 4 - actual_valence,
            Element::Nitrogen | Element::Phosphorus => {
                if actual_valence <= 3 {
                    3 - actual_valence
                }
                else if actual_valence <= 5 {
                    5 - actual_valence
                }
                else {
                    0
                }
            },
            Element::Oxygen => 2 - actual_valence,
            Element::Sulfur => {
                if actual_valence <= 2 {
                    2 - actual_valence
                }
                else if actual_valence <= 4 {
                    4 - actual_valence
                }
                else if actual_valence <= 6 {
                    6 - actual_valence
                }
                else {
                    0
                }
            },
            _ => 1 - actual_valence
        };
        for _ in 0..lack {
            let b = vec![sf.bonds.len()];
            sf.atoms[atom].bonds.push(sf.bonds.len());
            sf.bonds.push(Bond { a: atom, b: sf.atoms.len(), k: StructuralBond::Single });
            sf.atoms.push(AtomAndBondI { atom: Isotope::from(Element::Hydrogen), bonds: b });
        }
    }
    Ok(sf)
}