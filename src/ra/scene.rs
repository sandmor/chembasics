use std::mem;
use std::sync::Arc;

use rayon::prelude::*;
use image::{RgbImage, ImageBuffer};

use super::*;

const SHADOW_BIAS: PL = 1e-13;

#[derive(Debug)]
pub struct Scene {
    pub fov: f64,
    pub items: Vec<Item>,
    pub lights: Vec<Light>,
    pub width: usize,
    pub height: usize
}

impl Scene {
    fn trace<'a>(&'a self, ray: &Ray) -> Option<(PL, Color, &'a Item)> {
        let mut target = None;
        let mut d = PLMAX;
        for item in self.items.iter() {
            if let Some((od, c)) = item.intersect(&ray) {
                if od < d {
                    d = od;
                    target = Some((c, item));
                }
            }
        }
        if let Some((point_color, item)) = target {
            Some((d, point_color, item))
        }
        else {
            None
        }
    }

    fn get_color<'a>(&self, ray: &Ray, d: PL, obj_color: Color, item: &'a Item) -> Color {
        let hit_point = ray.orig + (ray.dir * d).as_point();
        let surface_normal = item.surface_normal(&hit_point);
        let reflection_ray = Ray {
            orig: hit_point + (surface_normal * SHADOW_BIAS),
            dir: ray.dir - (2.0 * ray.dir.dot(&surface_normal) * surface_normal)
        };
        let mut color = [0.; 3];
        for light in self.lights.iter() {
            let light_vector = (light.pos - hit_point).as_vector();
            let light_distance = light_vector.norm();
            let direction_to_light = light_vector.normalize();
            let shadow_ray = Ray {
                orig: hit_point + (direction_to_light * SHADOW_BIAS),
                dir: direction_to_light,
            };
            if let Some((distance, _, _)) = self.trace(&shadow_ray) {
                if distance < light_distance {
                    continue;
                }
            }
            let light_reflected = 1.0 / std::f32::consts::PI;
            let light_power = (light_distance / (4. * PLPI)) as f32 * light.intensity;
            let light_level = light_power * light_reflected;
            color[0] += obj_color[0] * light_level;
            color[1] += obj_color[1] * light_level;
            color[2] += obj_color[2] * light_level;
        }
        if color[0] < 0. {
            color[0] = 0.;
        }
        if color[1] < 0. {
            color[1] = 0.;
        }
        if color[2] < 0. {
            color[2] = 0.;
        }
        if color[0] > 1. {
            color[0] = 1.;
        }
        if color[1] > 1. {
            color[1] = 1.;
        }
        if color[2] > 1. {
            color[2] = 1.;
        }
        color
    }
}

pub fn draw(scene: Scene) {
    let mut img = vec![[0u8, 0u8, 0u8]; scene.width * scene.height];
    let aspect_ratio = (scene.width as f64) / (scene.height as f64);
    let scene = Arc::new(scene);
    img.par_iter_mut().enumerate().for_each(|(i, p)| {
        let y = i / scene.width;
        let x = i % scene.width;
        let ray = Ray::create_prime(aspect_ratio, x, y, &scene);
        if let Some((d, point_color, item)) = scene.trace(&ray) {
            let c = scene.get_color(&ray, d, point_color, item);
            *p = [(c[0] * 255.0) as u8, (c[1] * 255.0) as u8, (c[2] * 255.0) as u8];
        }
    });

    let mut plane_img = Vec::with_capacity(scene.width * scene.height * 3);

    for pixel in img {
        plane_img.push(pixel[0]);
        plane_img.push(pixel[1]);
        plane_img.push(pixel[2]);
    }

    let rimg: RgbImage = ImageBuffer::from_vec(scene.width as u32, scene.height as u32, plane_img).unwrap();
    rimg.save(concat!(env!("CARGO_MANIFEST_DIR"), "/result.png")).unwrap();
}