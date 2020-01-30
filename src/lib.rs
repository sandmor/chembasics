use std::mem;

use ptable::Element;

mod ion;
mod isotope;
mod basic_formulas;
mod structural;
pub mod consts;

pub use ion::*;
pub use isotope::*;
pub use structural::*;
pub use basic_formulas::*;

fn parse_element(string: &[u8]) -> (Option<Element>, &[u8]) {
    let mut end = 1;
    if string.len() >= 2 && string[1] >= b'a' && string[1] <= b'z' {
        end += 1;
    }
    (Element::from_symbol(unsafe { mem::transmute(&string[..end]) }), &string[end..])
}

fn parse_number(string: &[u8]) -> (usize, &[u8]) {
    let mut digits = Vec::new();
    let mut end = 0;
    while string.len() > end && (string[end] >= b'0' && string[end] <= b'9') {
        digits.push(string[end] - b'0');
        end += 1;
    }
    let mut result = 0;
    let mut m = 1;
    for d in digits.into_iter().rev() {
        result += d as usize * m;
        m *= 10;
    }
    (result, &string[end..])
}

pub trait AdvancedFormula {
    fn get_empirical_formula(&self) -> EmpiricalFormula;
}

pub trait BondClass {}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum QuantumBondKind {
    Sigma,
    Pi
}

impl BondClass for QuantumBondKind {}

#[cfg(test)]
mod tests {
    #[test]
    fn smiles_parse() {
        use crate::*;
        macro_rules! test {
            ($smiles:expr, $f:expr) => { assert_eq!(Compound::from_smiles(
                $smiles).unwrap().get_empirical_formula(), 
                EmpiricalFormula::from_string($f).unwrap()) };
        }
        test!("O", "H2O"); // Aqua
        test!("OO", "H2O2"); // Hydrogen peroxide
        test!("OS(=O)(=O)O", "H2SO4"); // Sulfuric acid
        test!("c1ccccc1", "H6C6"); // Benzene
        test!("CC(O)=O", "H4C2O2"); // Acetic acid
        test!("CC(C(=O)O)O", "C3H6O3"); // Lactic acid
        test!("CCOC(=O)C(=C)C#N", "C6H7NO2"); // Ethyl 2-cyanoacrylate
        test!("c1cc2c(c(c1)N)c(=O)[nH][nH]c2=O", "C8H7N3O2"); // Luminol
        test!("C1=CC(=CC=C1C(=O)O)C(=O)O.C1=CC(=CC=C1N)N", "C14H14N2O4"); // Poly(p-phenylene terephthalamide)
        test!("C[N+](C)(C)CCOP(=O)([O-])OCC(COC(=O)CCCCCCCCCCCOC(=O)CCCCC1CSSC1)OC(=O)CCCCCCCCCCCOC(=O)CCCCC2CSSC2", "C48H88NO12PS4"); // Dilipoyl lipid
        test!("[2H]", "H"); // Deuterium
    }
}
