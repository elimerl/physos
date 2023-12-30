use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;
use std::time::Duration;
use std::time::Instant;

use flo_canvas::*;
use flo_draw::*;
use flo_stream::*;

use futures::executor;
use futures::prelude::*;
use futures_timer::Delay;

use crate::math::Mat2;
use crate::math::Vec2;
use crate::World;

const SCALE: f32 = 50.0;

pub fn start_rendering(mut world: World) {
    with_2d_graphics(move || {
        // Create a window
        let (canvas, events) = create_drawing_window_with_events("Circle");
        let tick_stream = tick_stream(world.timestep);
        let events = events.map(|evt| WindowEvent::DrawEvent(evt));
        let mut events = stream::select(events, tick_stream);

        canvas.draw(|gc| {
            gc.clear_canvas(Color::Rgba(1.0, 1.0, 1.0, 1.0));
        });

        let mut canvas_size = Vec2::ONE;

        executor::block_on(async move {
            while let Some(event) = events.next().await {
                match event {
                    WindowEvent::Tick => {
                        let start = Instant::now();
                        world.step();
                        println!("tick took {:.4}s", (Instant::now() - start).as_secs_f32());

                        canvas.draw(|gc| {
                            gc.layer(LayerId(0));
                            gc.clear_layer();

                            for (idx, object) in &mut world.objects {
                                gc.new_path();
                                gc.identity_transform();
                                gc.canvas_height(canvas_size.y / SCALE);
                                gc.center_region(
                                    0.0,
                                    0.0,
                                    canvas_size.x / SCALE,
                                    canvas_size.y / SCALE,
                                );

                                // gc.transform(Transform2D::translate(
                                //     object.transform.position.x,
                                //     object.transform.position.y,
                                // ));
                                // gc.transform(Transform2D::rotate(object.transform.orientation));
                                match &object.shape {
                                    crate::Shape::Circle { radius } => {
                                        gc.circle(
                                            object.transform.position.x,
                                            object.transform.position.y,
                                            *radius,
                                        );
                                        gc.move_to(
                                            object.transform.position.x,
                                            object.transform.position.y,
                                        );
                                        gc.line_to(
                                            object.transform.position.x + *radius,
                                            object.transform.position.y,
                                        );
                                    }
                                    crate::Shape::Polygon(poly) => {
                                        let first = object.transform.transform(poly.vertices[0]);
                                        gc.move_to(first.x, first.y);
                                        for vertex in &poly.vertices {
                                            let p = object.transform.transform(*vertex);
                                            gc.line_to(p.x, p.y);
                                        }
                                        gc.close_path();
                                        let center = object.transform.transform(poly.center());
                                        gc.move_to(center.x, center.y);
                                        let support = object
                                            .transform
                                            .transform(poly.support(Vec2::new(1., 0.)));
                                        gc.line_to(support.x, support.y);

                                        gc.circle(
                                            object.transform.position.x,
                                            object.transform.position.y,
                                            poly.max_radius * 0.5,
                                        );
                                    }
                                }
                                let color = {
                                    let mut hasher = DefaultHasher::new();

                                    hasher.write_usize(idx.into_raw_parts().0);
                                    hasher.write_u64(idx.into_raw_parts().1);
                                    (hasher.finish() as f32 / 100.0).rem_euclid(360.)
                                };

                                gc.stroke_color(Color::Hsluv(color, 75.0, 65.0, 1.0));
                                gc.line_width_pixels(2.);
                                gc.stroke();
                                if object.transform.position.y < 0. {
                                    object.transform.position.y = canvas_size.y / SCALE + 10.0;
                                    object.velocity = Vec2::new(0., -5.0);
                                    object.angular_velocity = 0.;
                                }
                            }
                        });
                    }

                    WindowEvent::DrawEvent(DrawEvent::Closed) => {
                        break;
                    }

                    WindowEvent::DrawEvent(DrawEvent::Resize(x, y)) => {
                        canvas_size = Vec2::new(x as f32, y as f32);
                    }
                    _ => {}
                }
            }
        })
    });
}

fn tick_stream(timestep: f32) -> impl Send + Unpin + Stream<Item = WindowEvent> {
    generator_stream(move |yield_value| {
        async move {
            // Set up the clock
            let start_time = Instant::now();
            let mut last_time = Duration::from_millis(0);

            let max_ticks_per_call = 5;

            let tick_length = Duration::from_secs_f32(timestep);

            loop {
                let elapsed = start_time.elapsed() - last_time;

                let mut remaining = elapsed;
                let mut num_ticks = 0;
                while remaining >= tick_length {
                    if num_ticks < max_ticks_per_call {
                        yield_value(WindowEvent::Tick).await;
                        num_ticks += 1;
                    }

                    remaining -= tick_length;
                    last_time += tick_length;
                }

                // Wait for half a tick before generating more ticks
                let next_time = tick_length - remaining;
                let wait_time = Duration::min(tick_length / 2, next_time);

                Delay::new(wait_time).await;
            }
        }
        .boxed()
    })
}

enum WindowEvent {
    Tick,
    DrawEvent(DrawEvent),
}
