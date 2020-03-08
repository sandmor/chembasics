use std::num::{ParseIntError, ParseFloatError};
use std::str::Utf8Error;
use std::collections::BTreeMap;
use std::io::{self, Read};
use std::ops::Deref;

mod smiles;
mod cml;
mod mol;
mod nomenclature;
pub use mol::MolFile;

use crate::Point;
use crate::*;

include!(concat!(env!("OUT_DIR"), "/valences.rs"));

#[derive(Debug)]
pub enum ParserError {
    Syntax,
    IO(io::Error),
    Utf8(Utf8Error),
    UnexpectedEof,
    ParseInt(ParseIntError),
    ParseFloat(ParseFloatError),
}

impl From<xml::reader::Error> for ParserError {
    fn from(x: xml::reader::Error) -> ParserError {
        use xml::reader::ErrorKind;
        match x.kind().clone() {
            ErrorKind::Syntax(_) => ParserError::Syntax,
            ErrorKind::Io(e) => ParserError::IO(e),
            ErrorKind::Utf8(e) => ParserError::Utf8(e),
            ErrorKind::UnexpectedEof => ParserError::UnexpectedEof
        }
    }
}

impl From<ParseIntError> for ParserError {
    fn from(p: ParseIntError) -> ParserError {
        ParserError::ParseInt(p)
    }
}

impl From<ParseFloatError> for ParserError {
    fn from(p: ParseFloatError) -> ParserError {
        ParserError::ParseFloat(p)
    }
}

impl From<io::Error> for ParserError {
    fn from(e: io::Error) -> ParserError {
        ParserError::IO(e)
    }
}

// The order is important, that is used for optimize various things, DO NOT ALTER
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum StructuralBond {
    Aromatic,
    Single,
    Double,
    Triple,
}

impl BondClass for StructuralBond {}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Bond<B: BondClass> {
    pub a: usize,
    pub b: usize,
    pub k: B
}

impl<B: BondClass> Bond<B> {
    pub fn new(a: usize, b: usize, k: B) -> Bond<B> {
        Bond { a, b, k }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomAndBondI {
    pub atom: Isotope,
    pub bonds: Vec<usize>
}

impl Deref for AtomAndBondI {
    type Target = Isotope;
    
    fn deref(&self) -> &Isotope {
        &self.atom
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<AtomAndBondI>,
    pub bonds: Vec<Bond<StructuralBond>>,
    pub coords: Option<Vec<Point>>,
}

impl Molecule {
    pub fn from_smiles(string: &str) -> Result<Molecule, ()> {
        smiles::parse(string)
    }

    fn get_empirical_formula_optimize(&self, empirical: &mut BTreeMap<Element, usize>) {
        for a in self.atoms.iter() {
            *empirical.entry(*a.atom.get_element()).or_insert(0) += 1;
        }
    }

    pub fn atom_coords(&mut self) -> Option<&[Point]> {
        if let Some(ref r) = self.coords {
            Some(r)
        }
        else {
            None
        }
    }
}

impl BasicMolecule for Molecule {
    fn get_molecular_weight(&self) -> f32 {
        let mut weight = 0.0;
        for atom in self.atoms.iter() {
            weight += atom.atom.get_element().get_atomic_mass();
        }
        weight
    }
}

impl AdvancedFormula for Molecule {
    fn get_empirical_formula(&self) -> EmpiricalFormula {
        let mut empirical = BTreeMap::new();
        for a in self.atoms.iter() {
            *empirical.entry(a.atom).or_insert(0) += 1;
        }
        let mut res = Vec::with_capacity(empirical.len());
        for (k, v) in empirical {
            res.push((*k.get_element(), v));
        }
        EmpiricalFormula::new(res)
    }
}

/// A collection of multiple molecules
#[derive(Debug, Clone, PartialEq)]
pub struct Compound {
    molecules: Vec<Molecule>
}

impl Compound {
    pub fn from_smiles(smiles: &str) -> Result<Compound, ()> {
        let mut molecules = Vec::new();
        for molecule in smiles.split('.') {
            molecules.push(Molecule::from_smiles(molecule)?);
        }
        Ok(Compound { molecules })
    }

    pub fn from_cml<R: Read>(reader: R) -> Result<Compound, ParserError> {
        cml::parse(reader)
    }

    pub fn from_mol<R: Read>(reader: R) -> Result<Compound, ParserError> {
        Ok(MolFile::parse(reader)?.into_compound())
    }

    pub fn iter(&self) -> CompoundIterator {
        CompoundIterator { compound: &self, pos: 0 }
    }

    pub fn atoms_count(&self) -> usize {
        let mut total = 0;
        for molecule in self.molecules.iter() {
            total += molecule.atoms.len();
        }
        total
    }

    pub fn bonds_count(&self) -> usize {
        let mut total = 0;
        for molecule in self.molecules.iter() {
            total += molecule.bonds.len();
        }
        total
    }
}

impl BasicMolecule for Compound {
    fn get_molecular_weight(&self) -> f32 {
        let mut weight = 0.0;
        for g in self.molecules.iter() {
            weight += g.get_molecular_weight();
        }
        weight
    }
}

impl AdvancedFormula for Compound {
    fn get_empirical_formula(&self) -> EmpiricalFormula {
        let mut empirical = BTreeMap::new();
        for g in self.molecules.iter() {
            g.get_empirical_formula_optimize(&mut empirical);
        }
        let mut res = Vec::with_capacity(empirical.len());
        for (k, v) in empirical {
            res.push((k, v));
        }
        EmpiricalFormula::new(res)
    }
}

pub struct CompoundIterator<'a> {
    compound: &'a Compound,
    pos: usize
}

impl<'a> Iterator for CompoundIterator<'a> {
    type Item = &'a Molecule;
    fn next(&mut self) -> Option<Self::Item> {
        match self.compound.molecules.get(self.pos) {
            Some(m) => {
                self.pos += 1;
                Some(m)
            },
            None => {
                None
            }
        }
    }
}