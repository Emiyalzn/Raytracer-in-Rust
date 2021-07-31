use crate::bvh::*;
use crate::camera::*;
use crate::object::*;
use crate::vec3::*;
use rand::Rng;

fn random_double() -> f64 {
    rand::thread_rng().gen()
}

fn random_double_in(min: f64, max: f64) -> f64 {
    rand::thread_rng().gen_range(min, max)
}

pub fn init_scene(index: u32) -> (Arc<dyn Object>, Camera) {
    match index {
        1 => {
            let world_scene = two_spheres();
            let world = BvhNode::new_boxed(world_scene, 0.0, 0.0001);

            let look_from = Point3::new(13.0, 2.0, 3.0);
            let look_at = Point3::new(0.0, 0.0, 0.0);
            let vup = Vec3::new(0.0, 1.0, 0.0);
            let vfov = 20.0;
            let aspect_ratio = 3.0 / 2.0;
            let dist_to_focus = 10.0;
            let aperture = 0.0;
            let cam = Camera::new(
                vfov,
                aspect_ratio,
                look_from,
                look_at,
                vup,
                aperture,
                dist_to_focus,
            );

            (world, cam)
        }
        2 => {
            let world_scene = random_scene();
            let world = BvhNode::new_boxed(world_scene, 0.0, 0.0001);

            let look_from = Point3::new(13.0, 2.0, 3.0);
            let look_at = Point3::new(0.0, 0.0, 0.0);
            let vup = Vec3::new(0.0, 1.0, 0.0);
            let vfov = 20.0;
            let aspect_ratio = 3.0 / 2.0;
            let dist_to_focus = 10.0;
            let aperture = 0.1;
            let cam = Camera::new(
                vfov,
                aspect_ratio,
                look_from,
                look_at,
                vup,
                aperture,
                dist_to_focus,
            );

            (world, cam)
        }
        3 => {
            let world_scene = earth();
            let world = Arc::new(world_scene);

            let look_from = Point3::new(13.0, 2.0, 3.0);
            let look_at = Point3::new(0.0, 0.0, 0.0);
            let vup = Vec3::new(0.0, 1.0, 0.0);
            let vfov = 20.0;
            let aspect_ratio = 3.0 / 2.0;
            let dist_to_focus = 10.0;
            let aperture = 0.0;
            let cam = Camera::new(
                vfov,
                aspect_ratio,
                look_from,
                look_at,
                vup,
                aperture,
                dist_to_focus,
            );

            (world, cam)
        }
        _ => panic!("index out of bound"),
    }
}

pub fn earth() -> HittableList {
    let earth_texture = Arc::new(ImageTexture::new_by_pathstr("data/earthmap.jpg"));
    let earth_surface = Arc::new(Lambertian::new_arc(earth_texture));
    let globe = Arc::new(Sphere::new(Point3::new(0.0, 0.0, 0.0), 2.0, earth_surface));

    let mut world = HittableList::new();
    world.push(globe);
    world
}

pub fn two_spheres() -> HittableList {
    let mut objects = HittableList::new();
    let checker = Arc::new(CheckerTexture::new(
        &Color::new(0.2, 0.3, 0.1),
        &Color::new(0.9, 0.9, 0.9),
    ));
    let _checker = checker.clone();
    objects.push(Arc::new(Sphere::new(
        Point3::new(0.0, -10.0, 0.0),
        10.0,
        Arc::new(Lambertian::new_arc(checker)),
    )));
    objects.push(Arc::new(Sphere::new(
        Point3::new(0.0, 10.0, 0.0),
        10.0,
        Arc::new(Lambertian::new_arc(_checker)),
    )));
    objects
}

pub fn random_scene() -> HittableList {
    let mut world = HittableList::new();
    let checker_texture = Arc::new(CheckerTexture::new(
        &Color::new(0.2, 0.3, 0.1),
        &Color::new(0.9, 0.9, 0.9),
    ));
    world.push(Arc::new(Sphere::new(
        Point3::new(0.0, -1000.0, 0.0),
        1000.0,
        Arc::new(Lambertian::new_arc(checker_texture)),
    )));

    for a in -11..11 {
        for b in -11..11 {
            let choose_mat = random_double();
            let center = Point3::new(
                a as f64 + 0.9 * random_double(),
                0.2,
                b as f64 + 0.9 * random_double(),
            );

            if (center - Point3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                if choose_mat < 0.65 {
                    // diffuse
                    let albedo = Vec3::elemul(Color::random_unit(), Color::random_unit());
                    let sphere_material = Arc::new(Lambertian::new(&albedo));
                    world.push(Arc::new(Sphere::new(center, 0.2, sphere_material)));
                } else if choose_mat < 0.9 {
                    // metal
                    let albedo = Color::random(0.5, 1.0);
                    let fuzz = random_double_in(0.0, 0.5);
                    let sphere_material = Arc::new(Metal::new(&albedo, fuzz));
                    world.push(Arc::new(Sphere::new(center, 0.2, sphere_material)));
                } else {
                    // glass
                    let sphere_material = Arc::new(Dielectric::new(1.5));
                    world.push(Arc::new(Sphere::new(center, 0.2, sphere_material)));
                }
            }
        }
    }

    let material_1 = Arc::new(Dielectric::new(1.5));
    world.push(Arc::new(Sphere::new(
        Point3::new(0.0, 1.0, 0.0),
        1.0,
        material_1,
    )));

    let material_2 = Arc::new(Lambertian::new(&Color::new(0.4, 0.2, 0.1)));
    world.push(Arc::new(Sphere::new(
        Point3::new(-4.0, 1.0, 0.0),
        1.0,
        material_2,
    )));

    let material_3 = Arc::new(Metal::new(&Color::new(0.7, 0.6, 0.5), 0.0));
    world.push(Arc::new(Sphere::new(
        Point3::new(4.0, 1.0, 0.0),
        1.0,
        material_3,
    )));

    world
}
