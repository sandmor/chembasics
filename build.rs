use std::fs::File;
use std::path::Path;
use std::env;
use std::io::Write;

fn last_value(l: &[u8]) -> u8 {
    let mut i = l.len();
    for (p, n) in l.iter().enumerate() {
        if *n == 0 {
            i = p;
            break;
        }
    }
    if i == 0 {
        0
    }
    else {
        i -= 1;
        l[i]
    }
}

fn main() {
    let mut output = File::create(Path::new(&env::var("OUT_DIR").unwrap()).join("valences.rs")).unwrap();
    output.write(b"const VALENCES: [u8; 118] = [").unwrap();
    let mut first = true;
    for e in ptable::periodic_table() {
        if first {
            first = false;
        }
        else {
            output.write(b", ").unwrap();
        }
        let c = e.get_electronic_configuration();
        let mut v = 0;
        v += last_value(&c.s);
        v += last_value(&c.p);
        v += last_value(&c.d);
        v += last_value(&c.f);
        output.write(format!("{}", v).as_bytes()).unwrap();
    }
    output.write(b"];\n").unwrap();
}