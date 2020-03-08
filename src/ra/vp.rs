use std::ops::{Add, Sub, Mul, Div};

use crate::ra::scene::*;

pub type PL = f64;
pub const PLMAX: PL = std::f64::MAX;
pub const PLPI: PL = std::f64::consts::PI;

#[derive(Debug, Copy, Clone)]
pub struct Vector {
    pub x: PL,
    pub y: PL,
    pub z: PL,
}

impl Vector {
    pub fn new(x: PL, y: PL, z: PL) -> Vector {
        Vector { x, y, z }
    }

    pub fn as_point(self) -> Point {
        Point { x: self.x, y: self.y, z: self.z }
    }

    pub fn normalize(&self) -> Vector {
        let d = self.norm();
        Vector { x: self.x / d, y: self.y / d, z: self.z / d }
    }

    pub fn dot(self, other: &Vector) -> PL {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn norm(&self) -> PL {
        PL::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
    }
}

impl Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Vector {
        Vector { x: self.x + other.x, y: self.y + other.y, z: self.z + other.z }
    }
}

impl Sub<Vector> for Vector {
    type Output = Vector;

    fn sub(self, other: Vector) -> Vector {
        Vector { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
    }
}

impl Mul<PL> for Vector {
    type Output = Vector;
    fn mul(self, f: PL) -> Vector {
        Vector { x: self.x * f, y: self.y * f, z: self.z * f }
    }
}

impl Div<PL> for Vector {
    type Output = Vector;
    fn div(self, f: PL) -> Vector {
        Vector { x: self.x / f, y: self.y / f, z: self.z / f }
    }
}

impl Mul<Vector> for PL {
    type Output = Vector;
    fn mul(self, v: Vector) -> Vector {
        Vector { x: v.x * self, y: v.y * self, z: v.z * self }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point {
    pub x: PL,
    pub y: PL,
    pub z: PL,
}

impl Point {
    pub fn new(x: PL, y: PL, z: PL) -> Point {
        Point { x, y, z }
    }

    pub fn as_vector(self) -> Vector {
        Vector { x: self.x, y: self.y, z: self.z }
    }
}

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
    }
}

impl Add<Vector> for Point {
    type Output = Point;

    fn add(self, other: Vector) -> Point {
        Point { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
    }
}

impl Sub<Point> for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub orig: Point,
    pub dir: Vector,
}

impl Ray {
    pub fn create_prime(aspect_ratio: f64, x: usize, y: usize, scene: &Scene) -> Ray {
        let fov_adjustment = (scene.fov.to_radians() / 2.0).tan();
        let sensor_x = ((((x as f64 + 0.5) / scene.width as f64) * 2.0 - 1.0) * aspect_ratio) * fov_adjustment;
        let sensor_y = (1.0 - ((y as f64 + 0.5) / scene.height as f64) * 2.0) * fov_adjustment;

        Ray { orig: Point::new(0., 0., 0.), dir: Vector::new(sensor_x, sensor_y, -1.).normalize() }
    }
}