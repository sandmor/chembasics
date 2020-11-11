use std::io::{BufReader, BufRead, Read};

use fnv::FnvHashMap;
use ptable::Element;

use crate::Isotope;
use crate::ra::Point;
use super::{AtomAndBondI, Compound, Molecule, ParserError};

mod v2000;

#[derive(Debug)]
pub struct MolFile {
    name: String,
    comment: String,
    compound: Compound
}

impl MolFile {
    pub fn parse<R: Read>(reader: R) -> Result<MolFile, ParserError> {
        let mut reader = BufReader::new(reader).lines();
        macro_rules! get_or_eof {
            () => {{
                match reader.next() {
                    Some(line) => line?,
                    None => {
                        return Err(ParserError::UnexpectedEof);
                    }
                }
            }};
        }
        let name = get_or_eof!();
        let _ = get_or_eof!();
        let comment = get_or_eof!();
        let counts = get_or_eof!();
        if counts.len() < 39 {
            return Err(ParserError::Syntax);
        }
        let mut atoms_count = v2000::get_int_value(&counts[0..3])?;
        let mut bonds_count = v2000::get_int_value(&counts[3..6])?;
        let _atom_list_count = v2000::get_int_value(&counts[6..9])?;
        let mut chiral_flag = v2000::get_int_value(&counts[9..12])?;  // 0 or 1
        let _stext_entries_count = v2000::get_int_value(&counts[12..15])?;
        if chiral_flag > 1 {
            return Err(ParserError::Syntax);
        }
        let table_format = &counts[33..39];
        if table_format != "V2000" && table_format != "V2000" {
            return Err(ParserError::Syntax);
        }
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        if table_format == "V2000" {
            atoms.reserve(atoms_count as usize);
            bonds.reserve(bonds_count as usize);
            v2000::parse_atom_list(&mut reader, &mut atoms, atoms_count)?;
            v2000::parse_bond_list(&mut reader, &mut bonds, atoms_count as usize, bonds_count)?;
        }
        else {
            #[derive(Debug)]
            enum Block {
                Global,
                CTab,
                Atom,
                Bond,
            }
            let mut current_block = Block::Global;
            loop {
                let line = get_or_eof!();
                if line.len() < 6 {
                    return Err(ParserError::Syntax);
                }
                if !line.starts_with("M  ") {
                    return Err(ParserError::Syntax);
                }
                if &line[3..6] == "END" {
                    break;
                }
                if line.len() < 7 {
                    return Err(ParserError::Syntax);
                }
                if &line[3..7] != "V30 " {
                    return Err(ParserError::Syntax);
                }
                let mut line = (&line[7..]).split(' ');
                macro_rules! get_value {
                    () => {{
                        match line.next() {
                            Some(v) => v,
                            None => {
                                return Err(ParserError::Syntax);
                            }
                        }
                    }};
                }
                let head = get_value!();
                if head == "BEGIN" {
                    let key = get_value!();
                    match current_block {
                        Block::Global => {
                            if key == "CTAB" {
                                current_block = Block::CTab;
                            }
                        },
                        Block::CTab => {
                            if key == "ATOM" {
                                if atoms.len() != atoms_count as usize {
                                    return Err(ParserError::Syntax);
                                }
                                current_block = Block::Atom;
                            }
                            else if key == "BOND" {
                                current_block = Block::Bond;
                            }
                        }
                        _ => {}
                    }
                }
                else if head == "END" {
                    let key = get_value!();
                    match current_block {
                        Block::CTab => {
                            if key == "CTAB" {
                                current_block = Block::Global;
                            }
                        },
                        Block::Atom => {
                            if key == "ATOM" {
                                current_block = Block::CTab;
                            }
                        },
                        Block::Bond => {
                            if key == "BOND" {
                                current_block = Block::CTab;
                            }
                        },
                        _ => {}
                    }
                }
                if let Block::CTab = current_block {
                    if head == "COUNTS" {
                        atoms_count = get_value!().parse()?;
                        atoms.resize(atoms_count as usize, (Point::new(0.0, 0.0, 0.0), Element::Hydrogen));
                        bonds_count = get_value!().parse()?;
                        let _sgroups_count:u32 = get_value!().parse()?;
                        let _3d_constrants_count:u32 = get_value!().parse()?;
                        chiral_flag = get_value!().parse()?;
                        if chiral_flag > 1 {
                            return Err(ParserError::Syntax);
                        }
                    }
                }
                else if let Block::Atom = current_block {
                    let mut index:usize = get_value!().parse()?;
                    if index == 0 {
                        return Err(ParserError::Syntax);
                    }
                    index -= 1;
                    if index > atoms_count as usize {
                        return Err(ParserError::Syntax);
                    }
                    let el = match Element::from_symbol(get_value!()) {
                        Some(e) => e,
                        None => {
                            return Err(ParserError::Syntax);
                        }
                    };
                    let x:f64 = get_value!().parse()?;
                    let y:f64 = get_value!().parse()?;
                    let z:f64 = get_value!().parse()?;
                    atoms[index] = (Point::new(x, y, z), el);
                }
            }
        }
        let mut processed = FnvHashMap::with_capacity_and_hasher(atoms_count as usize, Default::default());
        let mut molecules:Vec<Molecule> = Vec::new();
        for bond in bonds.iter() {
            if let Some(i) = processed.get(&bond.a) {
                let (mol_id, atom_id):(usize, usize) = *i;
                if let None = processed.get(&bond.b) {
                    let bond_id = molecules[mol_id].bonds.len();
                    // Insert this bond on the molecule
                    molecules[mol_id].bonds.push(*bond);
                    // Update a atom bond list
                    molecules[mol_id].atoms[atom_id].bonds.push(bond_id);
                    // Get atom b id
                    let b_atom_id = molecules[mol_id].atoms.len();
                    // Register atom b
                    molecules[mol_id].atoms.push(AtomAndBondI { atom: Isotope::from(atoms[bond.b].1), bonds: vec![bond_id] });
                    // Register the atom b coords
                    molecules[mol_id].coords.as_mut().unwrap().push(atoms[bond.b].0);
                    processed.insert(bond.b, (mol_id, b_atom_id));
                }
            }
            else if let Some(i) = processed.get(&bond.b) {
                let (mol_id, atom_id):(usize, usize) = *i;
                let bond_id = molecules[mol_id].bonds.len();
                // Insert this bond on the molecule
                molecules[mol_id].bonds.push(*bond);
                // Update atom b bond list
                molecules[mol_id].atoms[atom_id].bonds.push(bond_id);
                // Get a atom id
                let a_atom_id = molecules[mol_id].atoms.len();
                // Register atom b
                molecules[mol_id].atoms.push(AtomAndBondI { atom: Isotope::from(atoms[bond.a].1), bonds: vec![bond_id] });
                // Register the atom a coords
                molecules[mol_id].coords.as_mut().unwrap().push(atoms[bond.a].0);
                processed.insert(bond.a, (mol_id, a_atom_id));
            }
            else {
                let mol_id = molecules.len();
                molecules.push(Molecule { atoms: vec![AtomAndBondI { 
                    atom: Isotope::from(atoms[bond.a].1), bonds: vec![0] },
                    AtomAndBondI { atom: Isotope::from(atoms[bond.b].1), bonds: vec![0] }], 
                    bonds: vec![*bond], coords: Some(vec![atoms[bond.a].0, atoms[bond.b].0]) });
                processed.insert(bond.a, (mol_id, 0));
                processed.insert(bond.b, (mol_id, 1));
            };
        }
        Ok(MolFile { comment, compound: Compound { molecules }, name })
    }

    pub fn comment(&self) -> &str {
        &self.comment
    }

    pub fn compound(&self) -> &Compound {
        &self.compound
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn set_comment(&mut self, comment: String) {
        self.comment = comment;
    }

    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    pub fn compound_mut(&mut self) -> &mut Compound {
        &mut self.compound
    }

    pub fn into_compound(self) -> Compound {
        self.compound
    }
}
