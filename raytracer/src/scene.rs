#![allow(dead_code)]
// You SHOULD remove above line in your code.

use crate::Vec3;
use rand::Rng;

fn random_double() -> f64 {
    rand::thread_rng().gen()
}

fn random_double_in(min: f64, max: f64) -> f64 {
    rand::thread_rng().gen_range(min, max)
}
