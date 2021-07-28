pub use crate::ray::Ray;
pub use crate::vec3::{Color, Point3, Vec3};
use rand::Rng;
pub use std::{sync::Arc, vec};

pub trait Object: Send + Sync {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

#[derive(Clone)]
pub struct HitRecord {
    pub p: Point3,
    pub normal: Vec3,
    pub t: f64,
    pub front_face: bool,
}

impl HitRecord {
    pub fn new(point: Point3, n: Vec3, tin: f64) -> Self {
        Self {
            p: Point3 {
                x: point.x,
                y: point.y,
                z: point.z,
            },
            normal: Vec3 {
                x: n.x,
                y: n.y,
                z: n.z,
            },
            t: tin,
            front_face: true,
        }
    }

    pub fn set_face_normal(&mut self, r: &Ray, outward_normal: &Vec3) {
        self.front_face = (r.dir * *outward_normal) < 0.0;
        self.normal = if self.front_face {
            *outward_normal
        } else {
            -*outward_normal
        };
    }
}

#[derive(Clone)]
pub struct Sphere {
    pub center: Point3,
    pub radius: f64,
}

impl Sphere {
    pub fn new(c: Point3, r: f64) -> Self {
        Self {
            center: Point3 {
                x: c.x,
                y: c.y,
                z: c.z,
            },
            radius: r,
        }
    }
}

impl Object for Sphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = r.orig - self.center;
        let a = r.dir.squared_length();
        let half_b = oc.dot(r.dir);
        let c = oc.squared_length() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;

        if discriminant > 0.0 {
            let root = discriminant.sqrt();
            let tmp = (-half_b - root) / a;
            if tmp < t_max && tmp > t_min {
                let n = (r.at(tmp) - self.center).unit();
                let outward_normal = (r.at(tmp) - self.center) / self.radius;
                let mut rec = HitRecord::new(r.at(tmp), n, tmp);
                rec.set_face_normal(r, &outward_normal);
                return Some(rec);
            }
            let tmp = (-half_b + root) / a;
            if tmp < t_max && tmp > t_min {
                let n = (r.at(tmp) - self.center).unit();
                let outward_normal = (r.at(tmp) - self.center) / self.radius;
                let mut rec = HitRecord::new(r.at(tmp), n, tmp);
                rec.set_face_normal(r, &outward_normal);
                return Some(rec);
            }
        }
        Option::None
    }
}

#[derive(Clone)]
pub struct HittableList {
    pub objects: Vec<Arc<dyn Object>>,
}

impl HittableList {
    pub fn new() -> Self {
        Self {
            objects: vec::Vec::new(),
        }
    }

    pub fn push(&mut self, ob: Arc<dyn Object>) {
        self.objects.push(ob);
    }
}

impl Object for HittableList {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest_so_far = t_max;
        let mut cur_rec = Option::<HitRecord>::None;

        for object in &self.objects {
            let this_rec = object.hit(ray, t_min, closest_so_far);
            if let Some(cur) = this_rec {
                closest_so_far = cur.t;
                cur_rec = Some(cur).clone();
            }
        }
        cur_rec
    }
}
