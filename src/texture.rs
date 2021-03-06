pub use crate::perlin::*;
use crate::vec3::*;
pub use image::*;
pub use std::path::*;
use std::sync::Arc;

pub trait Texture: Send + Sync {
    fn value(&self, u: f64, v: f64, p: &Point3) -> Color;
}

#[derive(Copy, Clone)]
pub struct SolidColor {
    pub color_value: Color,
}

impl SolidColor {
    pub fn new(c: &Color) -> Self {
        Self {
            color_value: Color {
                x: c.x,
                y: c.y,
                z: c.z,
            },
        }
    }

    pub fn new_rgb(red: f64, green: f64, blue: f64) -> Self {
        Self {
            color_value: Color {
                x: red,
                y: green,
                z: blue,
            },
        }
    }
}

impl Texture for SolidColor {
    fn value(&self, _u: f64, _v: f64, _p: &Point3) -> Color {
        self.color_value
    }
}

#[derive(Clone)]
pub struct CheckerTexture {
    odd: Arc<dyn Texture>,
    even: Arc<dyn Texture>,
}

impl CheckerTexture {
    pub fn new_arc(t0: Arc<dyn Texture>, t1: Arc<dyn Texture>) -> Self {
        Self { odd: t0, even: t1 }
    }

    pub fn new(c1: &Color, c2: &Color) -> Self {
        Self {
            odd: Arc::new(SolidColor::new(c1)),
            even: Arc::new(SolidColor::new(c2)),
        }
    }
}

impl Texture for CheckerTexture {
    fn value(&self, u: f64, v: f64, p: &Point3) -> Color {
        let sines = (10.0 * p.x).sin() * (10.0 * p.y).sin() * (10.0 * p.z).sin();
        if sines < 0.0 {
            self.odd.value(u, v, p)
        } else {
            self.even.value(u, v, p)
        }
    }
}

pub struct NoiseTexture {
    pub noise: Perlin,
}

impl NoiseTexture {
    pub fn new() -> Self {
        Self {
            noise: Perlin::new(),
        }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, _u: f64, _v: f64, p: &Vec3) -> Vec3 {
        return Vec3::ones() * self.noise.noise(p);
    }
}

impl Default for NoiseTexture {
    fn default() -> Self {
        Self::new()
    }
}

pub struct ImageTexture {
    pub data: ImageBuffer<image::Rgb<u8>, std::vec::Vec<u8>>,
}

fn clamp(x: f64, min: f64, max: f64) -> f64 {
    if x < min {
        return min;
    } else if x > max {
        return max;
    } else {
        return x;
    }
}

impl ImageTexture {
    pub fn new_by_pathstr(dir: &str) -> Self {
        return Self {
            data: image::open(&Path::new(dir)).unwrap().to_rgb(),
        };
    }

    pub fn width(&self) -> u32 {
        return self.data.width();
    }

    pub fn height(&self) -> u32 {
        return self.data.height();
    }
}

impl Texture for ImageTexture {
    fn value(&self, u: f64, v: f64, p: &Point3) -> Color {
        let _u = clamp(u, 0.0, 1.0);
        let _v = 1.0 - clamp(v, 0.0, 1.0);
        let mut i: u32 = (_u * self.width() as f64) as u32;
        let mut j: u32 = (_v * self.height() as f64) as u32;

        if i >= self.width() {
            i = self.width() - 1;
        }
        if j >= self.height() {
            j = self.height() - 1;
        }

        const COLOR_SCALE: f64 = 1.0 / 255.0;
        let pixel = self.data.get_pixel(i, j);
        let [red, green, blue] = pixel.0;
        return Color::new(
            red as f64 * COLOR_SCALE,
            green as f64 * COLOR_SCALE,
            blue as f64 * COLOR_SCALE,
        );
    }
}
