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
            let world_scene = random_scene_with_light();
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
        4 => {
            let world_scene = two_perlin_spheres();
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
        5 => {
            let world_scene = cornell_box();
            let world = BvhNode::new_boxed(world_scene, 0.0, 0.0001);

            let look_from = Point3::new(278.0, 278.0, -800.0);
            let look_at = Point3::new(278.0, 278.0, 0.0);
            let vfov = 40.0;
            let aspect_ratio = 1.0;
            let dist_to_focus = 10.0;
            let aperture = 0.0;
            let vup = Vec3::new(0.0, 1.0, 0.0);
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

pub fn cornell_box() -> HittableList {
    let red = Arc::new(Lambertian::new(&Color::new(0.65, 0.05, 0.05)));
    let white = Arc::new(Lambertian::new(&Color::new(0.73, 0.73, 0.73)));
    let green = Arc::new(Lambertian::new(&Color::new(0.12, 0.45, 0.15)));
    let light = Arc::new(DiffuseLight::new_color(&Color::new(15.0, 15.0, 15.0)));

    let mut world = HittableList::new();
    world.push(Arc::new(YzRect::new(0.0, 555.0, 0.0, 555.0, 555.0, green)));
    world.push(Arc::new(YzRect::new(0.0, 555.0, 0.0, 555.0, 0.0, red)));
    world.push(Arc::new(XzRect::new(
        213.0, 343.0, 227.0, 332.0, 554.0, light,
    )));
    world.push(Arc::new(XzRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        white.clone(),
    )));
    world.push(Arc::new(XzRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));
    world.push(Arc::new(XyRect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        white.clone(),
    )));

    let box1 = Arc::new(Box::new(
        &Point3::new(0.0, 0.0, 0.0),
        &Point3::new(165.0, 330.0, 165.0),
        white.clone(),
    ));
    let _box1 = Arc::new(RotateY::new(box1, 15.0));
    let _box1_ = Arc::new(Translate::new(_box1, &Vec3::new(265.0, 0.0, 295.0)));

    let box2 = Arc::new(Box::new(
        &Point3::new(0.0, 0.0, 0.0),
        &Point3::new(165.0, 165.0, 165.0),
        white.clone(),
    ));
    let _box2 = Arc::new(RotateY::new(box2, -18.0));
    let _box2_ = Arc::new(Translate::new(_box2, &Vec3::new(130.0, 0.0, 65.0)));

    world.push(_box1_);
    world.push(_box2_);

    world
}

pub fn two_perlin_spheres() -> HittableList {
    let mut world = HittableList::new();

    let pertext = Arc::new(NoiseTexture::new());
    let perlin_surface = Arc::new(Lambertian::new_arc(pertext));
    let _perlin_surface = perlin_surface.clone();

    world.push(Arc::new(Sphere::new(
        Point3::new(0.0, -1000.0, 0.0),
        1000.0,
        perlin_surface,
    )));
    world.push(Arc::new(Sphere::new(
        Point3::new(0.0, 2.0, 0.0),
        2.0,
        _perlin_surface,
    )));

    world
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

pub fn random_scene_with_light() -> HittableList {
    let mut world = HittableList::new();
    let checker = Arc::new(CheckerTexture::new(
        &Vec3::new(0.2, 0.3, 0.1),
        &Vec3::new(0.9, 0.9, 0.9),
    ));
    world.push(Arc::new(Sphere::new(
        Vec3::new(0.0, -1000.0, 0.0),
        -1000.0,
        Arc::new(Lambertian::new_arc(checker)),
    )));

    for a in -11..11 {
        for b in -11..11 {
            let choose_mat = random_double();
            let center = Vec3::new(
                a as f64 + 0.9 * random_double(),
                0.2,
                b as f64 + 0.9 * random_double(),
            );
            if (center - Vec3::new(4.0, 0.4, 0.0)).length() > 0.9 {
                let sphere_material: Arc<dyn Material>;

                if choose_mat < 0.8 {
                    // diffuse
                    let albedo = Vec3::elemul(Color::random_unit(), Color::random_unit());
                    sphere_material = Arc::new(Lambertian::new(&albedo));
                    world.push(Arc::new(Sphere::new(center, 0.2, sphere_material)));
                } else {
                    if choose_mat < 0.95 {
                        // metal
                        let albedo = Color::random(0.5, 1.0);
                        let sphere_material =
                            Arc::new(DiffuseLight::new(Arc::new(SolidColor::new(&albedo))));
                        world.push(Arc::new(Sphere::new(center, 0.2, sphere_material)));
                    } else {
                        // glass
                        sphere_material = Arc::new(Dielectric::new(1.5));
                        world.push(Arc::new(Sphere::new(center, 0.2, sphere_material)));
                    }
                }
            }
        }
    }
    let material1 = Arc::new(Dielectric::new(1.5));
    world.push(Arc::new(Sphere::new(
        Vec3::new(0.0, 1.0, 0.0),
        1.0,
        material1,
    )));

    let _checker = Arc::new(CheckerTexture::new(
        &Vec3::new(0.2, 0.3, 0.1),
        &Vec3::new(0.9, 0.9, 0.9),
    ));

    // light checker ball
    let material2 = Arc::new(DiffuseLight::new(_checker));
    world.push(Arc::new(Sphere::new(
        Vec3::new(4.0, 1.0, 0.0),
        1.0,
        material2,
    )));

    let material3 = Arc::new(Metal {
        albedo: Vec3::new(0.7, 0.6, 0.5),
        fuzz: 0.0,
    });
    world.push(Arc::new(Sphere::new(
        Vec3::new(-4.0, 1.0, 0.0),
        1.0,
        material3,
    )));

    return world;
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
