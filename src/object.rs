pub use crate::aabb::*;
pub use crate::ray::Ray;
pub use crate::texture::*;
pub use crate::vec3::{Color, Point3, Vec3};
use crate::{camera::degrees_to_radians, vec3::random_in_unit_sphere};
use rand::Rng;
pub use std::{sync::Arc, vec};

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
    pub u: f64,
    pub v: f64,
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
            u: 0.0,
            v: 0.0,
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

    pub fn set_uv(&mut self, res: (f64, f64)) {
        self.u = res.0;
        self.v = res.1;
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
                let res = get_sphere_uv(&rec.p);
                rec.set_uv(res);
                return Some(rec);
            }
            let tmp = (-half_b + root) / a;
            if tmp < t_max && tmp > t_min {
                let n = (r.at(tmp) - self.center).unit();
                let outward_normal = (r.at(tmp) - self.center) / self.radius;
                let mut rec = HitRecord::new(r.at(tmp), n, tmp, self.mat_ptr.clone());
                rec.set_face_normal(r, &outward_normal);
                let res = get_sphere_uv(&rec.p);
                rec.set_uv(res);
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

fn get_sphere_uv(p: &Vec3) -> (f64, f64) {
    let phi = p.z.atan2(p.x);
    let theta = p.y.asin();
    let u = 1.0 - (phi + std::f64::consts::PI) / (2.0 * std::f64::consts::PI);
    let v = (theta + std::f64::consts::PI / 2.0) / std::f64::consts::PI;
    (u, v)
}

pub struct XyRect {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub k: f64,
    pub mat_ptr: Arc<dyn Material>,
}

impl XyRect {
    pub fn new(x0: f64, x1: f64, y0: f64, y1: f64, k: f64, mat: Arc<dyn Material>) -> Self {
        Self {
            x0,
            x1,
            y0,
            y1,
            k,
            mat_ptr: mat,
        }
    }
}

impl Object for XyRect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - r.orig.z) / r.dir.z;
        if t < t_min || t > t_max {
            return None;
        }
        let x = r.orig.x + t * r.dir.x;
        let y = r.orig.y + t * r.dir.y;
        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None;
        }
        let res = (
            (x - self.x0) / (self.x1 - self.x0),
            (y - self.y0) / (self.y1 - self.y0),
        );
        let outward_normal = Vec3::new(0.0, 0.0, 1.0);
        let mut rec = HitRecord::new(r.at(t), outward_normal, t, self.mat_ptr.clone());
        rec.set_face_normal(r, &outward_normal);
        rec.set_uv(res);
        Some(rec)
    }
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let output_box = AABB::new(
            Vec3::new(self.x0, self.y0, self.k - 0.0001),
            Vec3::new(self.x1, self.y1, self.k + 0.0001),
        );
        Some(output_box)
    }
}

pub struct XzRect {
    pub x0: f64,
    pub x1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
    pub mat_ptr: Arc<dyn Material>,
}

impl XzRect {
    pub fn new(x0: f64, x1: f64, z0: f64, z1: f64, k: f64, mat: Arc<dyn Material>) -> Self {
        Self {
            x0,
            x1,
            z0,
            z1,
            k,
            mat_ptr: mat,
        }
    }
}

impl Object for XzRect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - r.orig.y) / r.dir.y;
        if t < t_min || t > t_max {
            return None;
        }
        let x = r.orig.x + t * r.dir.x;
        let z = r.orig.z + t * r.dir.z;
        if x < self.x0 || x > self.x1 || z < self.z0 || z > self.z1 {
            return None;
        }
        let res = (
            (x - self.x0) / (self.x1 - self.x0),
            (z - self.z0) / (self.z1 - self.z0),
        );
        let outward_normal = Vec3::new(0.0, 1.0, 0.0);
        let mut rec = HitRecord::new(r.at(t), outward_normal, t, self.mat_ptr.clone());
        rec.set_face_normal(r, &outward_normal);
        rec.set_uv(res);
        Some(rec)
    }
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let output_box = AABB::new(
            Vec3::new(self.x0, self.k - 0.0001, self.z0),
            Vec3::new(self.x1, self.k + 0.0001, self.z1),
        );
        Some(output_box)
    }
}

pub struct YzRect {
    pub y0: f64,
    pub y1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
    pub mat_ptr: Arc<dyn Material>,
}

impl YzRect {
    pub fn new(y0: f64, y1: f64, z0: f64, z1: f64, k: f64, mat: Arc<dyn Material>) -> Self {
        Self {
            y0,
            y1,
            z0,
            z1,
            k,
            mat_ptr: mat,
        }
    }
}

impl Object for YzRect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - r.orig.x) / r.dir.x;
        if t < t_min || t > t_max {
            return None;
        }
        let y = r.orig.y + t * r.dir.y;
        let z = r.orig.z + t * r.dir.z;
        if y < self.y0 || y > self.y1 || z < self.z0 || z > self.z1 {
            return None;
        }
        let res = (
            (y - self.y0) / (self.y1 - self.y0),
            (z - self.z0) / (self.z1 - self.z0),
        );
        let outward_normal = Vec3::new(1.0, 0.0, 0.0);
        let mut rec = HitRecord::new(r.at(t), outward_normal, t, self.mat_ptr.clone());
        rec.set_face_normal(r, &outward_normal);
        rec.set_uv(res);
        Some(rec)
    }
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let output_box = AABB::new(
            Vec3::new(self.k - 0.0001, self.y0, self.z0),
            Vec3::new(self.k + 0.0001, self.y1, self.z1),
        );
        Some(output_box)
    }
}

pub struct Box {
    box_min: Point3,
    box_max: Point3,
    sides: HittableList,
}

impl Box {
    pub fn new(p0: &Point3, p1: &Point3, mat: Arc<dyn Material>) -> Self {
        let mut objects = HittableList::new();
        objects.push(Arc::new(XyRect::new(
            p0.x,
            p1.x,
            p0.y,
            p1.y,
            p1.z,
            mat.clone(),
        )));
        objects.push(Arc::new(XyRect::new(
            p0.x,
            p1.x,
            p0.y,
            p1.y,
            p0.z,
            mat.clone(),
        )));
        objects.push(Arc::new(XzRect::new(
            p0.x,
            p1.x,
            p0.z,
            p1.z,
            p1.y,
            mat.clone(),
        )));
        objects.push(Arc::new(XzRect::new(
            p0.x,
            p1.x,
            p0.z,
            p1.z,
            p0.y,
            mat.clone(),
        )));
        objects.push(Arc::new(YzRect::new(
            p0.y,
            p1.y,
            p0.z,
            p1.z,
            p1.x,
            mat.clone(),
        )));
        objects.push(Arc::new(YzRect::new(
            p0.y,
            p1.y,
            p0.z,
            p1.z,
            p0.x,
            mat.clone(),
        )));

        Self {
            box_min: *p0,
            box_max: *p1,
            sides: objects,
        }
    }
}

impl Object for Box {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        return self.sides.hit(r, t_min, t_max);
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let output_box = AABB::new(self.box_min, self.box_max);
        Some(output_box)
    }
}

pub struct Translate {
    ptr: Arc<dyn Object>,
    offset: Vec3,
}

impl Translate {
    pub fn new(p: Arc<dyn Object>, displacement: &Vec3) -> Self {
        Self {
            ptr: p,
            offset: *displacement,
        }
    }
}

impl Object for Translate {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let moved_r = Ray::new(r.orig - self.offset, r.dir);
        match self.ptr.hit(&moved_r, t_min, t_max) {
            None => None,
            Some(mut rec) => {
                rec.p += self.offset;
                let normal = rec.normal;
                rec.set_face_normal(&moved_r, &normal);
                Some(rec)
            }
        }
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        match self.ptr.bounding_box(t0, t1) {
            None => None,
            Some(output_box) => {
                let new_box = AABB::new(
                    output_box.min_p + self.offset,
                    output_box.max_p + self.offset,
                );
                Some(new_box)
            }
        }
    }
}

pub struct RotateY {
    ptr: Arc<dyn Object>,
    sin_theta: f64,
    cos_theta: f64,
    hasbox: bool,
    bbox: AABB,
}

impl RotateY {
    pub fn new(p: Arc<dyn Object>, angle: f64) -> Self {
        let radians = degrees_to_radians(angle);
        let sine = radians.sin();
        let cosine = radians.cos();

        let mut min = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Point3::new(-f64::INFINITY, -f64::INFINITY, -f64::INFINITY);

        match p.bounding_box(0.0, 1.0) {
            None => Self {
                ptr: p,
                sin_theta: sine,
                cos_theta: cosine,
                hasbox: false,
                bbox: AABB::new(min, max),
            },
            Some(bound_box) => {
                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..2 {
                            let x =
                                i as f64 * bound_box.max_p.x + (1 - i) as f64 * bound_box.min_p.x;
                            let y =
                                j as f64 * bound_box.max_p.y + (1 - j) as f64 * bound_box.min_p.y;
                            let z =
                                k as f64 * bound_box.max_p.z + (1 - k) as f64 * bound_box.min_p.z;

                            let newx = cosine * x + sine * z;
                            let newz = -sine * x + cosine * z;

                            let tester = Vec3::new(newx, y, newz);

                            min.x = min.x.min(tester.x);
                            max.x = max.x.max(tester.x);
                            min.y = min.y.min(tester.y);
                            max.y = max.y.max(tester.y);
                            min.z = min.z.min(tester.z);
                            max.z = max.z.max(tester.z);
                        }
                    }
                }
                Self {
                    ptr: p,
                    sin_theta: sine,
                    cos_theta: cosine,
                    hasbox: true,
                    bbox: AABB::new(min, max),
                }
            }
        }
    }
}

impl Object for RotateY {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut orig = r.orig;
        let mut direction = r.dir;

        orig.x = self.cos_theta * r.orig.x - self.sin_theta * r.orig.z;
        orig.z = self.sin_theta * r.orig.x + self.cos_theta * r.orig.z;

        direction.x = self.cos_theta * r.dir.x - self.sin_theta * r.dir.z;
        direction.z = self.sin_theta * r.dir.x + self.cos_theta * r.dir.z;

        let rotated_r = Ray::new(orig, direction);

        match self.ptr.hit(&rotated_r, t_min, t_max) {
            None => None,
            Some(mut rec) => {
                let mut p = rec.p;
                let mut normal = rec.normal;

                p.x = self.cos_theta * rec.p.x + self.sin_theta * rec.p.z;
                p.z = -self.sin_theta * rec.p.x + self.cos_theta * rec.p.z;

                normal.x = self.cos_theta * rec.normal.x + self.sin_theta * rec.normal.z;
                normal.z = -self.sin_theta * rec.normal.x + self.cos_theta * rec.normal.z;

                rec.p = p;
                rec.set_face_normal(&rotated_r, &normal);
                Some(rec)
            }
        }
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        match self.hasbox {
            false => None,
            true => Some(self.bbox),
        }
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
    fn emitted(&self, u: f64, v: f64, p: &Point3) -> Color;
}

pub struct Lambertian {
    pub albedo: Arc<dyn Texture>,
}

impl Lambertian {
    pub fn new(a: &Color) -> Self {
        let tex = SolidColor::new(a);
        Self {
            albedo: Arc::new(tex),
        }
    }

    pub fn new_arc(a: Arc<dyn Texture>) -> Self {
        Self { albedo: a }
    }
}

impl Material for Lambertian {
    fn scatter(&self, r_in: &Ray, rec: &HitRecord) -> Option<(Color, Ray)> {
        let mut scatter_direction = rec.normal + crate::vec3::random_unit_vector();
        if scatter_direction.near_zero() {
            scatter_direction = rec.normal;
        }
        let scattered = Ray::new(rec.p, scatter_direction);
        let attenuation = self.albedo.value(rec.u, rec.v, &rec.p);
        Some((attenuation, scattered))
    }

    fn emitted(&self, u: f64, v: f64, p: &Point3) -> Color {
        Color::zero()
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

    fn emitted(&self, u: f64, v: f64, p: &Point3) -> Color {
        Color::zero()
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

    fn emitted(&self, u: f64, v: f64, p: &Point3) -> Color {
        Color::zero()
    }
}

pub struct DiffuseLight {
    emit: Arc<dyn Texture>,
}

impl DiffuseLight {
    pub fn new(a: Arc<dyn Texture>) -> Self {
        Self { emit: a }
    }

    pub fn new_color(c: &Color) -> Self {
        Self {
            emit: Arc::new(SolidColor::new(c)),
        }
    }
}

impl Material for DiffuseLight {
    fn scatter(&self, r_in: &Ray, rec: &HitRecord) -> Option<(Color, Ray)> {
        None
    }

    fn emitted(&self, u: f64, v: f64, p: &Point3) -> Color {
        self.emit.value(u, v, p)
    }
}
