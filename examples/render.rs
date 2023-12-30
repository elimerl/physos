use std::sync::Arc;

use nophysics::{math::Vec2, render::start_rendering, PhysicsObject, Polygon, Shape, World};

fn main() {
    let mut world = World::new(1. / 60., 3);
    let square = Arc::new(Polygon::rectangle(1., 1.));
    let floor = Arc::new(Polygon::rectangle(100., 1.));
    let wall = Arc::new(Polygon::rectangle(1., 100.));

    for i in 0..5 {
        world.objects.insert(
            PhysicsObject::new(Shape::polygon(Arc::clone(&square)), 1.)
                .with_position(Vec2::new(5.0, i as f32))
                .with_velocity(Vec2::new(0., 0.)),
        );
    }

    world.objects.insert(
        PhysicsObject::new(Shape::polygon(Arc::clone(&square)), 10.)
            .with_position(Vec2::new(10.0, 3.))
            .with_velocity(Vec2::new(-10., 0.)),
    );

    world.objects.insert(
        PhysicsObject::new(Shape::polygon(Arc::clone(&floor)), 1.)
            .with_frozen(true)
            .with_position(Vec2::new(10.0, 0.)),
    );
    world.objects.insert(
        PhysicsObject::new(Shape::polygon(Arc::clone(&wall)), 1.)
            .with_frozen(true)
            .with_position(Vec2::new(0.0, 0.)),
    );

    world.objects.insert(
        PhysicsObject::new(Shape::polygon(Arc::clone(&wall)), 1.)
            .with_frozen(true)
            .with_position(Vec2::new(20.0, 0.)),
    );

    start_rendering(world);
}
