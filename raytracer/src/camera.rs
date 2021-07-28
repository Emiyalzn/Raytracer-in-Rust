pub use crate::ray::Ray;
use crate::vec3::random_in_unit_disk;
pub use crate::vec3::Point3;
pub use crate::vec3::Vec3;

pub fn degrees_to_radians(degrees: f64) -> f64 {
    degrees * std::f64::consts::PI / 180.0
}

#[derive(Copy, Clone)]
pub struct Camera {
    pub origin: Point3,
    pub lower_left_corner: Point3,
    pub horizontal: Vec3,
    pub vertical: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub w: Vec3,
    pub lens_radius: f64,
}

impl Camera {
    pub fn new(
        vfov: f64,
        aspect_ratio: f64,
        look_from: Point3,
        look_at: Point3,
        vup: Vec3,
        aperture: f64,
        focus_dist: f64,
    ) -> Self {
        let theta = degrees_to_radians(vfov);
        let h = (theta / 2.0).tan();
        let viewport_height = 2.0 * h;
        let viewport_width = aspect_ratio * viewport_height;

        let _w = (look_from - look_at).unit();
        let _u = Vec3::cross(vup, _w).unit();
        let _v = Vec3::cross(_w, _u);

        Self {
            origin: look_from,
            horizontal: _u * viewport_width * focus_dist,
            vertical: _v * viewport_height * focus_dist,
            lower_left_corner: look_from
                - _u * viewport_width / 2.0 * focus_dist
                - _v * viewport_height / 2.0 * focus_dist
                - _w * focus_dist,
            u: _u,
            v: _v,
            w: _w,
            lens_radius: aperture / 2.0,
        }
    }

    pub fn get_ray(&self, s: f64, t: f64) -> Ray {
        let rd = crate::vec3::random_in_unit_disk() * self.lens_radius;
        let offset = self.u * rd.x + self.v * rd.y;
        Ray::new(
            self.origin + offset,
            self.lower_left_corner + self.horizontal * s + self.vertical * t - self.origin - offset,
        )
    }
}
