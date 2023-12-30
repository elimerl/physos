use std::{
    collections::{HashMap, HashSet},
    sync::Arc,
};

use generational_arena::{Arena, Index};
use math::{Mat2, Vec2};

pub mod math;
pub mod render;

pub const GRID_SIZE: f32 = 8.0;

/// A physics world.
#[derive(Clone, Debug)]
pub struct World {
    pub objects: Arena<PhysicsObject>,
    pub timestep: f32,
    pub gravity: Vec2,
    pub substeps: usize,

    changes: Vec<(Index, Vec2, Vec2)>,
}

#[derive(Clone, Copy, Debug, Default, Hash, PartialEq, Eq)]
struct GridPos(i32, i32);
impl GridPos {
    fn floor(pos: Vec2) -> GridPos {
        let pos = pos * (1. / GRID_SIZE);
        GridPos(pos.x.trunc() as i32, pos.y.trunc() as i32)
    }
}

impl World {
    pub fn new(timestep: f32, substeps: usize) -> Self {
        Self {
            objects: Arena::new(),
            gravity: Vec2::new(0.0, -10.0),
            timestep,
            changes: Vec::new(),
            substeps,
        }
    }

    fn step(&mut self) {
        let timestep = self.timestep / self.substeps as f32;
        for i in 0..self.substeps {
            for (idx, object) in &mut self.objects {
                if !object.frozen {
                    object.velocity += object.force * (1.0 / object.mass()) * timestep;
                    object.velocity += self.gravity * timestep;
                    object.transform.position += object.velocity * timestep;

                    object.angular_velocity += object.torque * object.inv_inertia() * timestep;
                    object.transform.orientation += object.angular_velocity * timestep;
                }
                object.force = Vec2::ZERO;
                object.torque = 0.;
            }

            for (idx, object) in &self.objects {
                if object.frozen {
                    continue;
                }
                for (other_idx, other_object) in &self.objects {
                    if other_idx == idx {
                        continue;
                    }

                    if !object.circle_overlap(other_object) {
                        continue;
                    }
                    if let Some((axis, overlap)) = object.overlap(other_object) {
                        let relative_velocity = object.velocity - other_object.velocity;
                        let v_j = -(1. + object.restitution) * relative_velocity.dot(axis);
                        let total_mass = object.inv_mass() + other_object.inv_mass();
                        let impulse = v_j / total_mass;
                        self.changes.push((
                            idx,
                            -axis * overlap * (object.inv_mass() / total_mass),
                            object.inv_mass() * impulse * axis,
                        ));
                    }
                }
            }
            for (idx, pos, vel) in &self.changes {
                let obj = self.objects.get_mut(*idx).unwrap();
                obj.transform.position += *pos;
                obj.velocity += *vel;
            }

            self.changes.clear();
        }
    }
}

#[derive(Clone, Debug)]
pub struct PhysicsObject {
    transform: Transform,
    velocity: Vec2,

    angular_velocity: f32,

    shape: Shape,

    density: f32,

    moment_of_inertia: f32,

    force: Vec2,
    torque: f32,

    frozen: bool,

    restitution: f32,
}

impl PhysicsObject {
    pub fn new(shape: Shape, density: f32) -> Self {
        let inertia = match &shape {
            Shape::Circle { radius } => 0.5 * shape.area() * density * radius * radius,
            Shape::Polygon(polygon) => {
                let mut inertia = 0.0;

                for face in polygon.faces() {
                    let a = face.0;
                    let b = face.1;
                    let mass_tri = density * 0.5 * a.cross(b).abs();
                    let inertia_tri =
                        mass_tri * (a.length_squared() + b.length_squared() + a.dot(b)) / 6.;
                    inertia += inertia_tri;
                }

                inertia
            }
        };
        Self {
            transform: Transform::identity(),
            shape,
            density,
            velocity: Vec2::ZERO,
            angular_velocity: 0.,
            moment_of_inertia: inertia,
            torque: 0.,
            restitution: 0.5,
            force: Vec2::ZERO,
            frozen: false,
        }
    }
    pub fn with_frozen(mut self, frozen: bool) -> PhysicsObject {
        self.frozen = frozen;
        self
    }
    pub fn with_position(mut self, position: Vec2) -> PhysicsObject {
        self.transform.position = position;
        self
    }
    pub fn with_velocity(mut self, velocity: Vec2) -> PhysicsObject {
        self.velocity = velocity;
        self
    }
    pub fn with_orientation(mut self, orientation: f32) -> PhysicsObject {
        self.transform.orientation = orientation;
        self
    }
    pub fn with_angular_velocity(mut self, angular_velocity: f32) -> PhysicsObject {
        self.angular_velocity = angular_velocity;
        self
    }
    pub fn with_torque(mut self, torque: f32) -> PhysicsObject {
        self.torque = torque;
        self
    }
    pub fn with_density(mut self, density: f32) -> PhysicsObject {
        self.density = density;
        self
    }
    pub fn with_shape(mut self, shape: Shape) -> PhysicsObject {
        self.shape = shape;
        self
    }

    pub fn inv_inertia(&self) -> f32 {
        1. / self.moment_of_inertia
    }

    pub fn inv_mass(&self) -> f32 {
        if self.frozen {
            0.
        } else {
            1. / self.mass()
        }
    }

    pub fn mass(&self) -> f32 {
        if self.frozen {
            f32::INFINITY
        } else {
            self.shape.area() * self.density
        }
    }

    pub fn add_force(&mut self, force: Vec2) {
        self.force += force;
    }
    pub fn add_torque(&mut self, torque: f32) {
        self.torque += torque;
    }

    pub fn add_force_at_pos(&mut self, force: Vec2, pos: Vec2) {
        self.add_force(force);
        self.add_torque(pos.cross(force));
    }

    pub fn add_impulse_at_pos(&mut self, impulse: Vec2, pos: Vec2) {
        self.velocity += impulse * self.inv_mass();
        self.angular_velocity += pos.cross(impulse) * self.inv_inertia();
    }

    fn overlap(&self, other: &PhysicsObject) -> Option<(Vec2, f32)> {
        match (&self.shape, &other.shape) {
            (
                Shape::Circle { radius },
                Shape::Circle {
                    radius: other_radius,
                },
            ) => todo!(),
            (Shape::Circle { radius }, Shape::Polygon(_)) => todo!(),
            (Shape::Polygon(_), Shape::Circle { radius }) => other.overlap(self),
            (Shape::Polygon(a), Shape::Polygon(b)) => {
                let mut depth = f32::INFINITY;
                let mut normal = Vec2::ZERO;
                for (v0, v1) in &a.faces {
                    let vec = self.transform.transform(*v1) - self.transform.transform(*v0);
                    let axis = Vec2::new(-vec.y, vec.x);
                    let proj_a = a.project(self.transform, axis);
                    let proj_b = b.project(other.transform, axis);

                    if proj_a.0 >= proj_b.1 || proj_b.0 >= proj_a.1 {
                        return None;
                    }

                    let axis_depth = f32::min(proj_b.1 - proj_a.0, proj_a.1 - proj_b.0);
                    if axis_depth < depth {
                        depth = axis_depth;
                        normal = axis;
                    }
                }
                for (v0, v1) in &b.faces {
                    let vec = self.transform.transform(*v1) - self.transform.transform(*v0);
                    let axis = Vec2::new(-vec.y, vec.x).normalized();
                    let proj_a = a.project(self.transform, axis);
                    let proj_b = b.project(other.transform, axis);

                    if proj_a.0 >= proj_b.1 || proj_b.0 >= proj_a.1 {
                        return None;
                    }

                    let axis_depth = f32::min(proj_b.1 - proj_a.0, proj_a.1 - proj_b.0);
                    if axis_depth < depth {
                        depth = axis_depth;
                        normal = axis;
                    }
                }

                if normal.dot(other.transform.position - self.transform.position) < 0. {
                    normal = -normal;
                }

                Some((normal, depth))
            }
        }
    }

    pub fn radius(&self) -> f32 {
        match &self.shape {
            Shape::Circle { radius } => *radius,
            Shape::Polygon(p) => p.max_radius,
        }
    }

    fn circle_overlap(&self, other_object: &PhysicsObject) -> bool {
        let r_a = self.radius();
        let r_b = other_object.radius();
        let dx = self.transform.position.x - other_object.transform.position.x;
        let dy = self.transform.position.y - other_object.transform.position.y;
        let r = r_a * r_b;
        dx * dx + dy * dy < r * r
    }
}

#[derive(Clone, Debug)]
pub enum Shape {
    Circle { radius: f32 },
    Polygon(Arc<Polygon>),
}

impl Shape {
    pub fn circle(radius: f32) -> Self {
        Self::Circle { radius }
    }
    pub fn polygon(polygon: Arc<Polygon>) -> Self {
        Self::Polygon(Arc::clone(&polygon))
    }

    pub fn area(&self) -> f32 {
        match &self {
            Shape::Circle { radius } => std::f32::consts::PI * radius * radius,
            Shape::Polygon(poly) => poly.area(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Polygon {
    vertices: Vec<Vec2>,
    faces: Vec<(Vec2, Vec2)>,
    center: Vec2,
    area: f32,
    face_normals: Vec<Vec2>,
    max_radius: f32,
}

impl Polygon {
    pub fn rectangle(width: f32, height: f32) -> Self {
        let left = -width / 2.;
        let right = left + width;
        let bottom = -height / 2.;
        let top = bottom + height;

        Self::new(vec![
            Vec2::new(left, top),
            Vec2::new(right, top),
            Vec2::new(right, bottom),
            Vec2::new(left, bottom),
        ])
    }

    pub fn new(mut vertices: Vec<Vec2>) -> Self {
        let center_avg = (
            vertices.iter().map(|v| v.x).sum::<f32>() / vertices.len() as f32,
            vertices.iter().map(|v| v.y).sum::<f32>() / vertices.len() as f32,
        );
        vertices.sort_by(|a, b| {
            let a = f32::atan2(a.x - center_avg.0, a.y - center_avg.1);
            let b = f32::atan2(b.x - center_avg.0, b.y - center_avg.1);
            b.total_cmp(&a) // clockwise
        });
        dbg!(&vertices);
        let faces = vertices
            .windows(2)
            .chain(std::iter::once(
                [vertices[vertices.len() - 1], vertices[0]].as_slice(),
            ))
            .map(|v| (v[0], v[1]))
            .collect::<Vec<_>>();

        let mut centroid = Vec2::ZERO;
        let polygon_area = faces
            .iter()
            .map(|face| 0.5 * face.0.cross(face.1).abs())
            .sum();
        for face in &faces {
            let face_area = 0.5 * face.0.cross(face.1).abs();
            centroid += (face.0 + face.1) * face_area / (3.0 * polygon_area);
        }

        for vertex in &mut vertices {
            *vertex -= centroid;
        }

        let face_normals = faces
            .iter()
            .map(|v| {
                let axis = v.1 - v.0;
                Vec2::new(-axis.y, axis.x)
            })
            .collect();

        let max_radius = vertices
            .iter()
            .map(|v| (centroid - *v).length())
            .max_by(|a, b| a.total_cmp(b))
            .unwrap()
            * 2.;

        Self {
            faces,
            face_normals,
            vertices,
            center: centroid,
            area: polygon_area,
            max_radius,
        }
    }

    pub fn support(&self, dir: Vec2) -> Vec2 {
        let mut best_proj = f32::NEG_INFINITY;
        let mut best_vert: Vec2 = Vec2::ZERO;

        for v in &self.vertices {
            let v = *v;
            let proj = v.dot(dir);
            if proj > best_proj {
                best_vert = v;
                best_proj = proj;
            }
        }

        best_vert
    }

    pub fn faces(&self) -> &[(Vec2, Vec2)] {
        &self.faces
    }
    pub fn face_normals(&self) -> &[Vec2] {
        &self.face_normals
    }

    pub fn center(&self) -> Vec2 {
        self.center
    }

    pub fn area(&self) -> f32 {
        self.area
    }

    pub fn project(&self, transform: Transform, axis: Vec2) -> (f32, f32) {
        let mut min = f32::INFINITY;
        let mut max = f32::NEG_INFINITY;

        for vertex in &self.vertices {
            let proj = transform.transform(*vertex).dot(axis);

            if proj < min {
                min = proj;
            } else if proj > max {
                max = proj;
            }
        }

        (min, max)
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Transform {
    position: Vec2,
    orientation: f32,
}
impl Transform {
    fn transform(&self, point: Vec2) -> Vec2 {
        (Mat2::rotation(self.orientation) * point) + self.position
    }
    fn identity() -> Self {
        Self::default()
    }
}

fn aabb_circle(c: Vec2, h: Vec2, p: Vec2, r: f32) -> bool {
    let v = p - c;
    let v = Vec2::new(v.x.abs(), v.y.abs());

    let u = v - h;
    let u = Vec2::new(u.x.max(0.), u.y.max(0.));
    u.length_squared() < r * r
}
