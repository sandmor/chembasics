use std::io::Read;
use std::hint::unreachable_unchecked;

use xml::{EventReader, reader::XmlEvent};
use fnv::FnvHashMap;
use ptable::Element;
use super::{AtomAndBondI, ParserError};
use crate::*;

#[derive(Debug)]
struct Data {
    atoms: FnvHashMap<String, (Element, Point)>,
    bond_index: FnvHashMap<String, Vec<usize>>,
    bonds: Vec<(String, String, StructuralBond)>
}

enum CMLRegion {
    AtomArray,
    BondArray,
    None
}

pub fn parse<R: Read>(reader: R) -> Result<Compound, ParserError> {
    let mut xml_reader = EventReader::new(reader);
    let mut in_molecule = false;
    let mut cml_region = CMLRegion::None;
    let mut atoms = FnvHashMap::default();
    let mut bonds = Vec::new();
    let mut bond_index = FnvHashMap::default();
    loop {
        let ev = xml_reader.next()?;
        match ev {
            XmlEvent::StartElement { name, attributes, namespace:_ } => {
                if name.namespace.is_some() {
                    // Ignore any unknown element
                    continue;
                }
                if !in_molecule {
                    if name.local_name == "molecule" {
                        in_molecule = true;
                    }
                    continue;
                }
                match cml_region {
                    CMLRegion::None => {
                        match &name.local_name[..] {
                            "atomArray" => {
                                cml_region = CMLRegion::AtomArray;
                            },
                            "bondArray" => {
                                cml_region = CMLRegion::BondArray;
                            },
                            _ => {}
                        }
                    },
                    CMLRegion::AtomArray => {
                        if name.local_name != "atom" {
                            continue;
                        }
                        let mut id = None;
                        let mut et = Element::Hydrogen;
                        let mut pos = Point::new(0., 0., 0.);
                        for attr in attributes.iter() {
                            if attr.name.namespace.is_some() {
                                continue;
                            }
                            match &attr.name.local_name[..] {
                                "id" => {
                                    id = Some(attr.value.to_owned());
                                },
                                "elementType" => {
                                    et = match Element::from_symbol(&attr.value) {
                                        Some(e) => e,
                                        None => {
                                            return Err(ParserError::Syntax);
                                        }
                                    };
                                },
                                "x3" => {
                                    pos.x = attr.value.parse()?;
                                },
                                "y3" => {
                                    pos.y = attr.value.parse()?;
                                },
                                "z3" => {
                                    pos.z = attr.value.parse()?;
                                },
                                _ => {}
                            }
                        }
                        let id = match id {
                            Some(id) => id,
                            None => {
                                return Err(ParserError::Syntax);
                            }
                        };
                        atoms.insert(id, (et, pos));
                    },
                    CMLRegion::BondArray => {
                        if name.local_name != "bond" {
                            continue;
                        }
                        let mut kind = StructuralBond::Single;
                        let mut a = String::new();
                        let mut b = String::new();
                        for attr in attributes.iter() {
                            if attr.name.namespace.is_some() {
                                continue;
                            }
                            match &attr.name.local_name[..] {
                                "atomRefs2" => {
                                    let mut ids = attr.value.split(' ');
                                    a = match ids.next() {
                                        Some(id) => id.to_owned(),
                                        None => {
                                            return Err(ParserError::Syntax);
                                        }
                                    };
                                    b = match ids.next() {
                                        Some(id) => id.to_owned(),
                                        None => {
                                            return Err(ParserError::Syntax);
                                        }
                                    };
                                    if let Some(_) = ids.next() {
                                        return Err(ParserError::Syntax);
                                    }
                                },
                                "order" => {
                                    match &attr.value[..] {
                                        "1" => {
                                            kind = StructuralBond::Single;
                                        },
                                        "2" => {
                                            kind = StructuralBond::Double;
                                        },
                                        "3" => {
                                            kind = StructuralBond::Triple;
                                        },
                                        _ => {
                                            return Err(ParserError::Syntax);
                                        }
                                    }
                                },
                                _ => {}
                            }
                        }
                        // If a not is empty then b tampoco
                        if a.is_empty() {
                            return Err(ParserError::Syntax);
                        }
                        bond_index.entry(a.clone()).or_insert(Vec::new()).push(bonds.len());
                        bond_index.entry(b.clone()).or_insert(Vec::new()).push(bonds.len());
                        bonds.push((a, b, kind));
                    },
                }
            },
            XmlEvent::EndElement { name } => {
                if name.namespace.is_some() {
                    // Ignore any unknown element
                    continue;
                }
                if in_molecule {
                    if name.local_name == "molecule" {
                        in_molecule = false;
                    }
                }
                match cml_region {
                    CMLRegion::AtomArray if name.local_name == "atomArray" => {
                        cml_region = CMLRegion::None;
                    },
                    CMLRegion::BondArray if name.local_name == "bondArray" => {
                        cml_region = CMLRegion::None;
                    },
                    _ => {}
                }
            },
            XmlEvent::EndDocument => {
                break;
            }
            _ => {}
        }
    }
    let mut molecules = Vec::new();
    let mut data = Data { atoms, bond_index, bonds };
    while !data.atoms.is_empty() {
        let (id, atom) = data.atoms.iter().next().unwrap();
        let id = id.to_owned();
        let atom = atom.to_owned();
        let mut molecule = Molecule { atoms: Vec::new(), bonds: Vec::new(), coords: Some(Vec::new()) };
        extract_group(&mut data, id, atom, &mut molecule);
        for (i, bond) in molecule.bonds.iter().enumerate() {
            molecule.atoms[bond.a].bonds.push(i);
            molecule.atoms[bond.b].bonds.push(i);
        }
        molecules.push(molecule);
    }
    Ok(Compound { molecules })
}

fn extract_group(data: &mut Data, id: String, atom: (Element, Point), group: &mut Molecule) {
    data.atoms.remove(&id);
    let my_molecule_id = group.atoms.len();
    group.atoms.push(AtomAndBondI { atom: Isotope::from(atom.0), bonds: Vec::new() });
    unsafe {
        match group.coords {
            Some(ref mut c) => {
                c.push(atom.1);
            },
            None => unreachable_unchecked()
        }
    }
    let bids = match data.bond_index.get(&id) {
        Some(bids) => bids.clone(),
        None => {
            return;
        }
    };
    for bid in bids {
        let bond = data.bonds[bid].clone();
        let pair;
        if bond.0 == id {
            pair = &bond.1;
        }
        else {
            pair = &bond.0;
        }
        if let Some(atom) = data.atoms.get(pair) {
            let atom = atom.clone();
            let next_atom_molecular_id = group.atoms.len();
            group.bonds.push(Bond { a: my_molecule_id, b: next_atom_molecular_id, k: bond.2 });
            extract_group(data, pair.to_owned(), atom, group);
        }
    }
}