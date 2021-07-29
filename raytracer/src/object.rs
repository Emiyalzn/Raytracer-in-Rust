pub use crate::ray::Ray;
pub use crate::aabb::*;
pub use crate::vec3::{Color, Point3, Vec3};
pub use std::{sync::Arc, vec};
use crate::vec3::random_in_unit_sphere;
use rand::Rng;

// trait的定义和使用方式？
pub trait Object: Send + Sync {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB>;
}

#[derive(Clone)]
pub struct HitRecord {
    pub p: Point3,
    pub normal: Vec3,
    pub mat_ptr: Arc<dyn Material>,
    pub t: f64,
    pub front_face: bool,
}

impl HitRecord {
    pub fn new(point: Point3, n: Vec3, tin: f64, m: Arc<dyn Material>) -> Self {
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
            mat_ptr: m,
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
    pub mat_ptr: Arc<dyn Material>,
}

impl Sphere {
    pub fn new(c: Point3, r: f64, m: Arc<dyn Material>) -> Self {
        Self {
            center: Point3 {
                x: c.x,
                y: c.y,
                z: c.z,
            },
            radius: r,
            mat_ptr: m,
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
                let mut rec = HitRecord::new(r.at(tmp), n, tmp, self.mat_ptr.clone());
                rec.set_face_normal(r, &outward_normal);
                return Some(rec);
            }
            let tmp = (-half_b + root) / a;
            if tmp < t_max && tmp > t_min {
                let n = (r.at(tmp) - self.center).unit();
                let outward_normal = (r.at(tmp) - self.center) / self.radius;
                let mut rec = HitRecord::new(r.at(tmp), n, tmp, self.mat_ptr.clone());
                rec.set_face_normal(r, &outward_normal);
                return Some(rec);
            }
        }
        Option::None
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        let output_box = AABB::new(
            self.center - Vec3::ones() * self.radius,
            self.center + Vec3::ones() * self.radius,
        );
        Some(output_box)
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

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        if self.objects.is_empty() {
            return None;
        }

        let mut output_box = AABB::new(Point3::zero(), Point3::zero());
        let mut first_box = true;

        for object in &self.objects {
            let res = object.bounding_box(t0, t1);
            match res {
                Some(tmp_box) => {
                    output_box = if first_box {
                        tmp_box
                    } else {
                        surrounding_box(output_box, tmp_box)
                    };
                    first_box = false;
                }
                None => {
                    return None;
                }
            }
        }
        Some(output_box)
    }
}

pub trait Material: Send + Sync {
    fn scatter(&self, r_in: &Ray, rec: &HitRecord) -> Option<(Color, Ray)>;
}

pub struct Lambertian {
    pub albedo: Color,
}

impl Lambertian {
    pub fn new(a: &Color) -> Self {
        Self {
            albedo: Color {
                x: a.x,
                y: a.y,
                z: a.z,
            },
        }
    }
}

impl Material for Lambertian {
    fn scatter(&self, r_in: &Ray, rec: &HitRecord) -> Option<(Color, Ray)> {
        let mut scatter_direction = rec.normal + crate::vec3::random_unit_vector();
        if scatter_direction.near_zero() {
            scatter_direction = rec.normal;
        }
        let scattered = Ray::new(rec.p, scatter_direction);
        let attenuation = self.albedo;
        Some((attenuation, scattered))
    }
}

pub struct Metal {
    pub albedo: Color,
    pub fuzz: f64,
}

impl Metal {
    pub fn new(a: &Color, f: f64) -> Self {
        Self {
            albedo: Color {
                x: a.x,
                y: a.y,
                z: a.z,
            },
            fuzz: if f < 1.0 { f } else { 1.0 },
        }
    }
}

impl Material for Metal {
    fn scatter(&self, r_in: &Ray, rec: &HitRecord) -> Option<(Color, Ray)> {
        let reflected = crate::vec3::reflect(&r_in.dir.unit(), &rec.normal);
        let scattered = Ray::new(rec.p, reflected + random_in_unit_sphere() * self.fuzz);
        let attenuation = self.albedo;
        if scattered.dir * rec.normal > 0.0 {
            Some((attenuation, scattered))
        } else {
            None
        }
    }
}

pub struct Dielectric {
    ir: f64,
}

impl Dielectric {
    pub fn new(ri: f64) -> Self {
        Self { ir: ri }
    }
}

pub fn reflectance(cosine: f64, ref_idx: f64) -> f64 {
    let mut r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * (1.0 - cosine).powf(5.0);
}

impl Material for Dielectric {
    fn scatter(&self, r_in: &Ray, rec: &HitRecord) -> Option<(Color, Ray)> {
        let attenuation = Color::ones();
        let refraction_ratio = if rec.front_face {
            1.0 / self.ir
        } else {
            self.ir
        };
        let unit_direction = r_in.dir.unit();
        let cos = -unit_direction * rec.normal;
        let cos_theta = if cos < 1.0 { cos } else { 1.0 };
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        let cannot_refract = refraction_ratio * sin_theta > 1.0;
        let mut direction = Vec3::zero();

        if cannot_refract || reflectance(cos_theta, refraction_ratio) > rand::thread_rng().gen() {
            direction = crate::vec3::reflect(&unit_direction, &rec.normal);
        } else {
            direction = crate::vec3::refract(&unit_direction, &rec.normal, refraction_ratio);
        }

        let scattered = Ray::new(rec.p, direction);
        Some((attenuation, scattered))
    }
}
