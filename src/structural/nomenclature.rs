use ptable::Element;
use fnv::FnvHashMap;

#[derive(Debug, Clone)]
struct ChainItem {
    id: usize,
    branches: Vec<Vec<ChainItem>>
}

fn name_branch(branch: Vec<ChainItem>) -> Option<String> {
    let mut name;
    let length = branch.len();
    let decene = length / 10;
    let unit = length % 10;
    if decene == 0 {
        name = match unit {
            1 => "met",
            2 => "et",
            3 => "prop",
            4 => "but",
            5 => "pent",
            6 => "hexa",
            7 => "hept",
            8 => "oct",
            9 => "non",
            _ => unreachable!("The unit is not a unit")
        }.to_owned();
    }
    else {
        name = match unit {
            0 if decene > 1 => "ei",
            0 => "",
            1 => "un",
            2 => "do",
            3 => "tri",
            4 => "tetra",
            5 => "penta",
            6 => "hexa",
            7 => "hept",
            8 => "oct",
            9 => "non",
            _ => unreachable!("The unit is not a unit")
        }.to_owned();
        if decene == 1 {
            name.push_str("dec");
        }
        else if decene == 2 {
            name.push_str("cos");
        }
        else {
            unimplemented!("Alkane too large");
        }
    }
    Some(name)
}

impl super::Molecule {
    fn seek_largest_chain(&self, last: Option<usize>, id: usize, among_chain: bool) -> Vec<ChainItem> {
        let mut chain;
        let mut branches = Vec::new();
        for bond in self.atoms[id].bonds.iter() {
            let bond = self.bonds[*bond];
            let pair;
            if bond.a == id {
                pair = bond.b;
            }
            else {
                pair = bond.a;
            }
            if *self.atoms[pair].atom.get_element() == Element::Hydrogen {
                // Ignore hydrogen atoms
                continue;
            }
            if Some(pair) == last {
                // Skip previous one
                continue;
            }
            let path = self.seek_largest_chain(Some(id), pair, false);
            branches.push(path);
        }
        if among_chain {
            chain = Vec::with_capacity(branches[0].len()+1+branches[1].len());
            chain.extend_from_slice(&branches[0]);
            chain.push(ChainItem { id, branches: Vec::new() });
            chain.extend_from_slice(&branches[1]);
        }
        else if branches.is_empty() {
            chain = vec![ChainItem { id, branches: Vec::new() }];
        }
        else {
            branches.sort_unstable_by_key(|b| b.len());
            chain = vec![ChainItem { id, branches: branches[1..].to_vec() }];
            chain.extend_from_slice(&branches[0]);
        }
        chain
    }

    pub fn get_name(&self) -> Option<String> {
        let mut name = String::new();
        for (id, atom) in self.atoms.iter().enumerate() {
            if *atom.get_element() == Element::Carbon {
                let mut bonds_count = 0;
                for bond in self.atoms[id].bonds.iter() {
                    let bond = self.bonds[*bond];
                    let pair;
                    if bond.a == id {
                        pair = bond.b;
                    }
                    else {
                        pair = bond.a;
                    }
                    if *self.atoms[pair].get_element() == Element::Hydrogen {
                        // Ignore hydrogen atoms
                        continue;
                    }
                    bonds_count += 1;
                }
                if bonds_count <= 2 {
                    if bonds_count == 0 {
                        name = "methane".to_owned();
                    }
                    else {
                        let chain;
                        if bonds_count == 2 {
                            chain = self.seek_largest_chain(None, id, true);
                        }
                        else {
                            chain = self.seek_largest_chain(None, id, false);
                        }
                        println!("{:?}", chain.len());
                        name = name_branch(chain)?;
                        if name.ends_with("et") {
                            name.push_str("hane");
                        }
                        else if name.ends_with('a') {
                            name.push_str("ne");
                        }
                        else {
                            name.push_str("ane");
                        }
                    }
                    break;
                }
            }
        }
        println!("{:?}", name);
        Some(name)
    }
}