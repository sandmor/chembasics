use crate::ra::Point;
use ptable::Element;
use crate::structural::*;
use std::io::{BufReader, Lines, Read};

pub fn parse_atom_list<R: Read>(lines: &mut Lines<BufReader<R>>, atoms_buffer: &mut Vec<(Point, Element)>, 
    atoms_count: u32) -> Result<(), ParserError> {
    for _ in 0..atoms_count {
        let line = match lines.next() {
            Some(line) => line?,
            None => {
                return Err(ParserError::UnexpectedEof);
            }
        };
        if line.len() < 69 {
            return Err(ParserError::Syntax);
        }
        let x = get_float_value(&line[0..10])?;
        let y = get_float_value(&line[10..20])?;
        let z = get_float_value(&line[20..30])?;
        let el = match Element::from_symbol((&line[31..34]).trim()) {
            Some(e) => e,
            None => {
                return Err(ParserError::Syntax);
            }
        };
        atoms_buffer.push((Point::new(x, y, z), el));
    }
    Ok(())
}

pub fn parse_bond_list<R: Read>(lines: &mut Lines<BufReader<R>>, bonds_buffer: &mut Vec<Bond<StructuralBond>>, 
    atoms_count: usize, bonds_count: u32) -> Result<(), ParserError> {
    for _ in 0..bonds_count {
        let line = match lines.next() {
            Some(line) => line?,
            None => {
                return Err(ParserError::UnexpectedEof);
            }
        };
        let mut line = line.split(' ');
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
        let mut a:usize = get_value!().parse()?;
        if a == 0 {
            return Err(ParserError::Syntax);
        }
        a -= 1;
        if a >= atoms_count {
            return Err(ParserError::Syntax);
        }
        let mut b:usize = get_value!().parse()?;
        if b == 0 {
            return Err(ParserError::Syntax);
        }
        b -= 1;
        if b >= atoms_count {
            return Err(ParserError::Syntax);
        }
        let k:usize = get_value!().parse()?;
        if k == 0 || k > 3 {
            return Err(ParserError::Syntax);
        }
        bonds_buffer.push(Bond::new(a, b, unsafe { mem::transmute::<u8, StructuralBond>(k as u8) }));
    }
    Ok(())
}

pub fn get_float_value(text: &str) -> Result<f64, ParserError> {
    if text.find(|c:char| !c.is_whitespace()) == None  {
        return Ok(0.0);
    }
    Ok(text.trim().parse()?)
}

pub fn get_int_value(text: &str) -> Result<u32, ParserError> {
    if text.find(|c:char| !c.is_whitespace()) == None  {
        return Ok(0);
    }
    Ok(text.trim().parse()?)
}