mod smiles;

use std::collections::{BTreeMap, BTreeSet};

use crate::*;

//include!(concat!(env!("OUT_DIR"), "/valences.rs"));   // This adds a VALENCE constant, which is an arrangement of valence electrons number of all periodic table atoms

// The order is important, that is used for optimize various things, DO NOT ALTER
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum StructuralBondKind {
    Aromatic,
    Single,
    Double,
    Triple,
}

impl BondClass for StructuralBondKind {}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Bond<B: BondClass> {
    pub a: usize,
    pub b: usize,
    pub k: B
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Molecule {
    pub atoms: Vec<(Ion, Vec<usize>)>,
    pub bonds: Vec<Bond<StructuralBondKind>>,
}

impl Molecule {
    pub fn from_smiles(string: &str) -> Result<Molecule, ()> {
        let mut sf = Molecule { atoms: Vec::new(), bonds: Vec::new() };
        let string = string.as_bytes();
        let mut misc = smiles::SMILESMisc { automatic_hydrogens_targets: Vec::new(), 
            labels: BTreeMap::new(), aromatic_ions: BTreeSet::new() };
        smiles::parse_smiles_group(string, &mut sf, &mut misc, None, &mut smiles::AromaticDetectionData { init: false, 
            last_one_was_double: false }, &mut Vec::new())?;
        for atom in misc.automatic_hydrogens_targets {
            let mut actual_valence = 0;
            let mut aromatic_flip_flop = true;
            for bond in sf.atoms[atom].1.iter() {
                if sf.bonds[*bond].k == StructuralBondKind::Aromatic {
                    actual_valence += 1;
                    if aromatic_flip_flop {
                        actual_valence += 1;
                    }
                    aromatic_flip_flop = !aromatic_flip_flop;
                }
                actual_valence += sf.bonds[*bond].k as isize;
            }
            let lack = match sf.atoms[atom].0.get_element() {
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
                sf.atoms[atom].1.push(sf.bonds.len());
                sf.bonds.push(Bond { a: atom, b: sf.atoms.len(), k: StructuralBondKind::Single });
                sf.atoms.push((Ion::new(Element::Hydrogen, 0), b));
            }
        }
        Ok(sf)
    }

    fn get_empirical_formula_optimize(&self, empirical: &mut BTreeMap<Element, usize>) {
        for a in self.atoms.iter() {
            *empirical.entry(*a.0.get_element()).or_insert(0) += 1;
        }
    }
}

impl BasicMolecule for Molecule {
    fn get_molecular_weight(&self) -> f32 {
        let mut weight = 0.0;
        for (e, _) in self.atoms.iter() {
            weight += e.get_element().get_atomic_mass();
        }
        weight
    }
}

impl AdvancedFormula for Molecule {
    fn get_empirical_formula(&self) -> EmpiricalFormula {
        let mut empirical = BTreeMap::new();
        for a in self.atoms.iter() {
            *empirical.entry(a.0).or_insert(0) += 1;
        }
        let mut res = Vec::with_capacity(empirical.len());
        for (k, v) in empirical {
            res.push((*k.get_element(), v));
        }
        EmpiricalFormula::new(res)
    }
}

/// A collection of multiple molecules
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Compound {
    groups: Vec<Molecule>
}

impl Compound {
    pub fn from_smiles(smiles: &str) -> Result<Compound, ()> {
        let mut groups = Vec::new();
        for group in smiles.split('.') {
            groups.push(Molecule::from_smiles(group)?);
        }
        Ok(Compound { groups })
    }
}

impl BasicMolecule for Compound {
    fn get_molecular_weight(&self) -> f32 {
        let mut weight = 0.0;
        for g in self.groups.iter() {
            weight += g.get_molecular_weight();
        }
        weight
    }
}

impl AdvancedFormula for Compound {
    fn get_empirical_formula(&self) -> EmpiricalFormula {
        let mut empirical = BTreeMap::new();
        for g in self.groups.iter() {
            g.get_empirical_formula_optimize(&mut empirical);
        }
        let mut res = Vec::with_capacity(empirical.len());
        for (k, v) in empirical {
            res.push((k, v));
        }
        EmpiricalFormula::new(res)
    }
}