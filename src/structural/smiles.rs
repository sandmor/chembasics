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
pub fn parse_smiles_group<'a>(mut string: &'a [u8], sf: &mut Molecule, misc: &mut SMILESMisc, 
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
                        b':' => StructuralBondKind::Aromatic,
                        b'-' => StructuralBondKind::Single,
                        b'=' => StructuralBondKind::Double,
                        b'#' => StructuralBondKind::Triple,
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
                        k = StructuralBondKind::Aromatic;
                    }
                    else {
                        k = StructuralBondKind::Single;
                    }
                }
                current_bounds.push(sf.bonds.len());
                sf.atoms[adj_atom].1.push(sf.bonds.len());
                if aromatic_detection_data.init {
                    // Default believe there is no aromaticity
                    aromatic_detection_data.init = false;
                    if k == StructuralBondKind::Single || k == StructuralBondKind::Double {
                        if ((k == StructuralBondKind::Single) == aromatic_detection_data.last_one_was_double) {
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
                    if k == StructuralBondKind::Single || k == StructuralBondKind::Double {
                        aromatic_detection_data.init = true;
                        aromatic_detection_data.last_one_was_double = k == StructuralBondKind::Double;
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
                                StructuralBondKind::Aromatic
                            }
                            else {
                                StructuralBondKind::Single
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
                                    sf.bonds[i].k = StructuralBondKind::Aromatic;
                                }
                                k = StructuralBondKind::Aromatic;
                                bonds_in_the_probably_aromatic_ring.clear();
                                aromatic_detection_data.init = false;
                            }
                        }
                    }
                    current_bounds.push(sf.bonds.len());
                    sf.atoms[b].1.push(sf.bonds.len());
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
            sf.atoms.push((Isotope::new(Ion::new($atom, $charge), $isotopic_spec), current_bounds));
        }};
    }
    if let Some(_) = adj_atom {
        if !string.is_empty() {
            match string[0] {
                b'-' => {
                    waiting_bond = Some(StructuralBondKind::Single);
                    string = &string[1..];
                },
                b'=' => {
                    waiting_bond = Some(StructuralBondKind::Double);
                    string = &string[1..];
                },
                b'#' => {
                    waiting_bond = Some(StructuralBondKind::Triple);
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
                        waiting_bond = Some(StructuralBondKind::Aromatic);
                        continue;
                    },
                    b'-' => {
                        waiting_bond = Some(StructuralBondKind::Single);
                        continue;
                    },
                    b'=' => {
                        waiting_bond = Some(StructuralBondKind::Double);
                        continue;
                    },
                    b'#' => {
                        waiting_bond = Some(StructuralBondKind::Triple);
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