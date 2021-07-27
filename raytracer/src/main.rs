#![allow(clippy::float_cmp)]
#![feature(box_syntax)]

// mod, pub use, use的区别？

mod material;
mod ray;
mod scene;
mod vec3;

use image::{ImageBuffer, Rgb, RgbImage};
use indicatif::ProgressBar;
pub use ray::Ray;
use rusttype::Font;
use scene::example_scene;
use std::sync::mpsc::channel;
use std::sync::Arc;
use threadpool::ThreadPool;
pub use vec3::{Color, Point3, Vec3};

const AUTHOR: &str = "Emiyalzn";

pub struct World {
    pub height: u32,
}

impl World {
    pub fn color(&self, _: u32, y: u32) -> u8 {
        (y * 256 / self.height) as u8
    }
}

fn get_text() -> String {
    // GITHUB_SHA is the associated commit ID
    // only available on GitHub Action
    let github_sha = option_env!("GITHUB_SHA")
        .map(|x| "@".to_owned() + &x[0..6])
        .unwrap_or_default();
    format!("{}{}", AUTHOR, github_sha)
}

fn is_ci() -> bool {
    option_env!("CI").unwrap_or_default() == "true"
}

fn render_text(image: &mut RgbImage, msg: &str) {
    let font_file = if is_ci() {
        "EncodeSans-Regular.ttf"
    } else {
        "/System/Library/Fonts/Helvetica.ttc"
    };
    let font_path = std::env::current_dir().unwrap().join(font_file);
    let data = std::fs::read(&font_path).unwrap();
    let font: Font = Font::try_from_vec(data).unwrap_or_else(|| {
        panic!(format!(
            "error constructing a Font from data at {:?}",
            font_path
        ));
    });

    imageproc::drawing::draw_text_mut(
        image,
        Rgb([255, 255, 255]),
        10,
        10,
        rusttype::Scale::uniform(24.0),
        &font,
        msg,
    );
}

fn main() {
    // get environment variable CI, which is true for GitHub Action
    let is_ci = is_ci();

    // jobs: split image into how many parts
    // workers: maximum allowed concurrent running threads
    let (n_jobs, n_workers): (usize, usize) = if is_ci { (32, 2) } else { (16, 2) };

    println!(
        "CI: {}, using {} jobs and {} workers",
        is_ci, n_jobs, n_workers
    );

    // Image
    let aspect_ratio = 16.0 / 9.0;
    let image_width = 400;
    let image_height = (image_width as f64 / aspect_ratio) as u32;

    // Camera
    let viewport_height = 2.0;
    let viewport_width = aspect_ratio * viewport_height;
    let focal_length = 1.0;

    let origin = Color::new(0.0, 0.0, 0.0);
    let horizontal = Vec3::new(viewport_width, 0.0, 0.0);
    let vertical = Vec3::new(0.0, viewport_height, 0.0);
    let lower_left_corner =
        origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, focal_length);

    // create a channel to send objects between threads
    let (tx, rx) = channel();
    let pool = ThreadPool::new(n_workers);

    let bar = ProgressBar::new(n_jobs as u64);

    // use Arc to pass one instance of World to multiple threads
    let world = Arc::new(example_scene());

    for i in 0..n_jobs {
        let tx = tx.clone();
        let world_ptr = world.clone();
        pool.execute(move || {
            // here, we render some of the rows of image in one thread
            let row_begin = image_height as usize * i / n_jobs;
            let row_end = image_height as usize * (i + 1) / n_jobs;
            let render_height = row_end - row_begin;
            let mut img: RgbImage = ImageBuffer::new(image_width, render_height as u32);
            for x in 0..image_width {
                // img_y is the row in partial rendered image
                // y is real position in final image
                for (img_y, y) in (row_begin..row_end).enumerate() {
                    let y = y as u32;
                    let u: f64 = (x as f64) / (image_width - 1) as f64;
                    let v: f64 = (y as f64) / (image_height - 1) as f64;
                    let ray = Ray::new(
                        origin,
                        lower_left_corner + horizontal * u + vertical * v - origin,
                    );
                    let cur_color = ray_color(&ray);
                    write_color(cur_color, &mut img, x, img_y as u32);
                    // let pixel = img.get_pixel_mut(x, img_y as u32);
                    // let color = world_ptr.color(x, y);
                    // *pixel = Rgb([color, color, color]);
                }
            }
            // send row range and rendered image to main thread
            tx.send((row_begin..row_end, img))
                .expect("failed to send result");
        });
    }

    let mut result: RgbImage = ImageBuffer::new(image_width, image_height);

    for (rows, data) in rx.iter().take(n_jobs) {
        // idx is the corrsponding row in partial-rendered image
        for (idx, row) in rows.enumerate() {
            for col in 0..image_width {
                let row = row as u32;
                let idx = idx as u32;
                *result.get_pixel_mut(col, row) = *data.get_pixel(col, idx);
            }
        }
        bar.inc(1);
    }
    bar.finish();

    // render commit ID and author name on image
    let msg = get_text();
    println!("Extra Info: {}", msg);

    render_text(&mut result, msg.as_str());

    result.save("output/test.png").unwrap();
}

fn hit_sphere(center: Point3, radius: f64, r: &Ray) -> bool {
    let oc = r.orig - center;
    let a = r.dir.dot(r.dir);
    let b = oc.dot(r.dir) * 2.0;
    let c = oc.dot(oc) - radius * radius;
    let discriminant = b * b - 4.0 * a * c;
    discriminant > 0.0
}

fn ray_color(r: &Ray) -> Color {
    if hit_sphere(Point3::new(0.0, 0.0, -1.0), 0.5, r) {
        return Color::new(1.0, 0.0, 0.0);
    }
    let unit_direction = r.dir.unit();
    let t = 0.5 * (unit_direction.y + 1.0);
    Color {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    } * (1.0 - t)
        + Color {
            x: 0.5,
            y: 0.7,
            z: 1.0,
        } * t
}

fn within(min: f64, max: f64, value: f64) -> f64 {
    if value > max {
        return max;
    }
    if value < min {
        return min;
    }
    value
}

fn write_color(color: Color, img: &mut RgbImage, x: u32, y: u32) {
    let pixel = img.get_pixel_mut(x, y);
    let colorx = color.x.sqrt();
    let colory = color.y.sqrt();
    let colorz = color.z.sqrt();
    *pixel = image::Rgb([
        (within(0.0, 0.999, colorx) * 256.0) as u8,
        (within(0.0, 0.999, colory) * 256.0) as u8,
        (within(0.0, 0.999, colorz) * 256.0) as u8,
    ]);
}
