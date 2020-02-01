mod vp;
mod scene;

pub use vp::*;
pub use scene::*;

pub type Color = [f32; 3];

#[derive(Debug)]
pub enum ItemKind {
    Atom(f64, Color)
}

#[derive(Debug)]
pub struct Item {
    pub pos: Point,
    pub kind: ItemKind,
}

impl Item {
    pub fn intersect(&self, ray: &Ray) -> Option<(PL, Color)> {
        match self.kind {
            ItemKind::Atom(r, c) => {
                // Line segment between the ray origin and the atom center
                let l = (self.pos - ray.orig).as_vector();
                let adj_len = l.dot(&ray.dir);

                let d2 = l.dot(&l) - (adj_len * adj_len);  // This is equivalent to: l.length() * l.length() - adj_len * adj_len
                let radius2 = r * r;
                if d2 > radius2 {
                    return None;
                }
                let thc = (radius2 - d2).sqrt();
                let t0 = adj_len - thc;
                let t1 = adj_len + thc;
                if t0 < 0.0 && t1 < 0.0 {
                    return None;
                }
                Some((PL::min(t0, t1), c))
            }
        }
    }

    fn surface_normal(&self, hit_point: &Point) -> Vector {
        match self.kind {
            ItemKind::Atom(_, _) => (*hit_point - self.pos).as_vector().normalize(),
        }
    }
}

#[derive(Debug)]
pub struct Light {
    pub pos: Point,
    pub intensity: f32
}