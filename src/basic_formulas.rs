use std::iter::Iterator;
use std::ops::{Index, IndexMut};
use std::fmt::{self, Display};
use std::collections::BTreeMap;

use ptable::Element;

use crate::*;

pub trait BasicMolecule {
    fn get_molecular_weight(&self) -> f32;
}

#[derive(Debug, Clone, Ord, PartialOrd, Hash)]
pub struct EmpiricalFormula(Vec<(Element, usize)>);

impl EmpiricalFormula {
    pub fn new(formula: Vec<(Element, usize)>) -> EmpiricalFormula {
        EmpiricalFormula(formula)
    }

    pub fn from_string(string: &str) -> Result<EmpiricalFormula, ()> {
        let mut string = string.as_bytes();
        let mut elements = Vec::new();
        if string.is_empty() {
            return Err(());
        }
        while !string.is_empty() {
            if string[0] < b'A' || string[0] > b'Z' {
                return Err(());
            }
            let (e, s) = parse_element(string);
            let e = match e {
                Some(e) => e,
                None => {
                    return Err(());
                },
            };
            let mut count = 1;
            string = s;
            if s.is_empty() {
                elements.push((e, count));
                break;
            }
            if string[0] >= b'0' && string[0] <= b'9' {
                let (c, s) = parse_number(string);
                count = c;
                elements.push((e, count));
                string = s;
                if string.is_empty() {
                    break;
                }
            }
            else {
                elements.push((e, count));
            }
        }
        Ok(EmpiricalFormula(elements))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> EmpiricalFormulaIterator {
        EmpiricalFormulaIterator { items: self, count: 0 }
    }

    pub fn get(&self, i: usize) -> Option<&(Element, usize)> {
        self.0.get(i)
    }

    pub fn get_mut(&mut self, i: usize) -> Option<&mut (Element, usize)> {
        self.0.get_mut(i)
    }

    pub unsafe fn get_unchecked(&self, i: usize) -> &(Element, usize) {
        self.0.get_unchecked(i)
    }

    pub unsafe fn get_unchecked_mut(&mut self, i: usize) -> &mut (Element, usize) {
        self.0.get_unchecked_mut(i)
    }
}

impl BasicMolecule for EmpiricalFormula {
    fn get_molecular_weight(&self) -> f32 {
        let mut weight = 0.0;
        for (i, c) in self.0.iter() {
            weight += i.get_atomic_mass() * *c as f32;
        }
        weight
    }
}

impl Index<usize> for EmpiricalFormula {
    type Output = (Element, usize);

    fn index(&self, i: usize) -> &Self::Output {
        &self.0[i]
    }
}

impl IndexMut<usize> for EmpiricalFormula {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.0[i]
    }
}

impl Display for EmpiricalFormula {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for e in self.0.iter() {
            write!(f, "{}", e.0.get_symbol())?;
            if e.1 != 1 {
                write!(f, "{}", e.1)?;
            }
        }
        Ok(())
    }
}

pub struct EmpiricalFormulaIterator<'a> {
    items: &'a EmpiricalFormula,
    count: usize
}

impl<'a> Iterator for EmpiricalFormulaIterator<'a> {
    type Item = &'a (Element, usize);
    fn next(&mut self)  -> Option<Self::Item> {
        if self.count == self.items.len() {
            return None;
        }
        let item = &self.items[self.count];
        self.count += 1;
        Some(item)
    }
}

impl PartialEq for EmpiricalFormula {
    fn eq(&self, other: &Self) -> bool {
        let mut check = vec![false; self.0.len()];
        for (e, n) in other.0.iter() {
            let mut exist = false;
            for (i, (_e, _n)) in self.0.iter().enumerate() {
                if e == _e {
                    if n != _n {
                        return false;
                    }
                    check[i] = true;
                    exist = true;
                    break;
                }
            }
            if !exist {
                return false;
            }
        }
        for c in check {
            if c == false {
                return false;
            }
        }
        true
    }
}

impl Eq for EmpiricalFormula {}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ElementOrGroup {
    Element(Element),
    Group(MolecularFormula),
}

fn parse_group(mut string: &[u8]) -> Result<(MolecularFormula, &[u8]), ()> {
    if string.is_empty() {
        return Err(());
    }
    let mut result = Vec::new();
    loop {
        if string[0] == b')' {
            string = &string[1..];
            break;
        }
        let e = if string[0] == b'(' {
            let (g, s) = parse_group(&string[1..])?;
            string = s;
            ElementOrGroup::Group(g)
        }
        else {
            if string[0] < b'A' || string[0] > b'Z' {
                return Err(());
            }
            let (e, s) = parse_element(string);
            let e = match e {
                Some(e) => e,
                None => {
                    return Err(());
                },
            };
            string = s;
            ElementOrGroup::Element(e)
        };
        let mut count = 1;
        if string.is_empty() {
            return Err(());
        }
        if string[0] >= b'0' && string[0] <= b'9' {
            let (c, s) = parse_number(string);
            count = c;
            result.push((e, count));
            string = s;
            if string.is_empty() {
                return Err(());
            }
        }
        else {
            result.push((e, count));
        }
    }
    Ok((MolecularFormula(result), string))
}

#[derive(Debug, Clone, Ord, PartialOrd, Hash)]
pub struct MolecularFormula(Vec<(ElementOrGroup, usize)>);

impl MolecularFormula {
    pub fn new(formula: Vec<(ElementOrGroup, usize)>) -> MolecularFormula {
        MolecularFormula(formula)
    }

    pub fn from_string(string: &str) -> Result<MolecularFormula, ()> {
        let mut string = string.as_bytes();
        let mut result = Vec::new();
        if string.is_empty() {
            return Err(());
        }
        while !string.is_empty() {
            let e = if string[0] == b'(' {
                let (g, s) = parse_group(&string[1..])?;
                string = s;
                ElementOrGroup::Group(g)
            }
            else {
                if string[0] < b'A' || string[0] > b'Z' {
                    return Err(());
                }
                let (e, s) = parse_element(string);
                let e = match e {
                    Some(e) => e,
                    None => {
                        return Err(());
                    },
                };
                string = s;
                ElementOrGroup::Element(e)
            };
            let mut count = 1;
            if string.is_empty() {
                result.push((e, count));
                break;
            }
            if string[0] >= b'0' && string[0] <= b'9' {
                let (c, s) = parse_number(string);
                count = c;
                result.push((e, count));
                string = s;
                if string.is_empty() {
                    break;
                }
            }
            else {
                result.push((e, count));
            }
        }
        Ok(MolecularFormula(result))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> MolecularFormulaIterator {
        MolecularFormulaIterator { items: self, count: 0 }
    }

    pub fn get(&self, i: usize) -> Option<&(ElementOrGroup, usize)> {
        self.0.get(i)
    }

    pub fn get_mut(&mut self, i: usize) -> Option<&mut (ElementOrGroup, usize)> {
        self.0.get_mut(i)
    }

    pub unsafe fn get_unchecked(&self, i: usize) -> &(ElementOrGroup, usize) {
        self.0.get_unchecked(i)
    }

    pub unsafe fn get_unchecked_mut(&mut self, i: usize) -> &mut (ElementOrGroup, usize) {
        self.0.get_unchecked_mut(i)
    }

    fn get_empirical_formula_inner(&self, r: &mut BTreeMap<Element, usize>) {
        for item in self.0.iter() {
            match item.0 {
                ElementOrGroup::Element(e) => {
                    *r.entry(e).or_insert(0) += 1;
                },
                ElementOrGroup::Group(ref g) => {
                    g.get_empirical_formula_inner(r);
                },
            }
        }
    }

    pub fn get_empirical_formula(&self) -> EmpiricalFormula {
        let mut r = BTreeMap::new();
        self.get_empirical_formula_inner(&mut r);
        let mut res = Vec::new();
        for item in r {
            res.push(item);
        }
        EmpiricalFormula::new(res)
    }

    fn get_empirical_formula_optimize(&self, empirical: &mut BTreeMap<Element, usize>) {
        for a in self.0.iter() {
            match a.0 {
                ElementOrGroup::Element(e) =>  {
                    *empirical.entry(e).or_insert(0) += a.1;
                },
                ElementOrGroup::Group(ref g) => {
                    g.get_empirical_formula_optimize(empirical);
                }
            }
        }
    }
}

impl BasicMolecule for MolecularFormula {
    fn get_molecular_weight(&self) -> f32 {
        let mut weight = 0.0;
        for (i, c) in self.0.iter() {
            match i {
                ElementOrGroup::Element(e) => {
                    weight += e.get_atomic_mass() * *c as f32;
                },
                ElementOrGroup::Group(g) => {
                    weight += g.get_molecular_weight();
                }
            }
        }
        weight
    }
}

impl Index<usize> for MolecularFormula {
    type Output = (ElementOrGroup, usize);

    fn index(&self, i: usize) -> &Self::Output {
        &self.0[i]
    }
}

impl IndexMut<usize> for MolecularFormula {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.0[i]
    }
}

pub struct MolecularFormulaIterator<'a> {
    items: &'a MolecularFormula,
    count: usize
}

impl<'a> Iterator for MolecularFormulaIterator<'a> {
    type Item = &'a (ElementOrGroup, usize);
    fn next(&mut self)  -> Option<Self::Item> {
        if self.count == self.items.len() {
            return None;
        }
        let item = &self.items[self.count];
        self.count += 1;
        Some(item)
    }
}

impl PartialEq for MolecularFormula {
    fn eq(&self, other: &Self) -> bool {
        let mut check = vec![false; self.0.len()];
        for (e, n) in other.0.iter() {
            let mut exist = false;
            for (i, (_e, _n)) in self.0.iter().enumerate() {
                if e == _e {
                    if n != _n {
                        return false;
                    }
                    check[i] = true;
                    exist = true;
                    break;
                }
            }
            if !exist {
                return false;
            }
        }
        for c in check {
            if c == false {
                return false;
            }
        }
        true
    }
}

impl Eq for MolecularFormula {}

impl AdvancedFormula for MolecularFormula {
    fn get_empirical_formula(&self) -> EmpiricalFormula {
        let mut empirical = BTreeMap::new();
        self.get_empirical_formula_optimize(&mut empirical);
        let mut res = Vec::with_capacity(empirical.len());
        for (k, v) in empirical {
            res.push((k, v));
        }
        EmpiricalFormula::new(res)
    }
}
