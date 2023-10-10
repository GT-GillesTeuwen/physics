use std::ptr::null;

use sdl2::gfx::primitives::DrawRenderer;
use sdl2::keyboard::Scancode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;
use sdl2::render::{Canvas, TextureQuery};
use sdl2::ttf::Font;
use sdl2::video::Window;
use sdl2::EventPump;

mod events;
mod physics;
mod render;

#[derive(Clone, Copy)]
pub struct Vector2 {
    x: f32,
    y: f32,
}
#[derive(Clone, Copy)]
pub struct RigidBody {
    mass: f32,
    velocity: Vector2,
    acceleration: Vector2,
}

#[derive(Clone)]
pub struct History {
    trail: Vec<Vector2>,
}
#[derive(Clone)]
pub struct SimulationObject {
    shape: Circle,
    selected: bool,
    position: Vector2,
    body: RigidBody,
    material: Material,
    history: History,
}

#[derive(Clone, Copy)]
pub struct Material {
    restitution: f32,
    colour: Color,
}

#[derive(Clone, Copy)]
pub struct Circle {
    radius: f32,
}

fn main() {
    let ttf_context = sdl2::ttf::init().unwrap();

    let font_path = "C:/Windows/Fonts/AGENCYR.TTF"; // Replace with your actual font path
    let font = ttf_context.load_font(font_path, 24).unwrap();

    let mut is_paused = true;
    let mut single_step = false;

    let dampening = 0.7;
    let screen_width = 1000;
    let screen_height = 800;
    let mut selected_index = 0;
    let force_strength = 8000.0; // Now in units per second
    let r = 40;

    let gravitational_constant = 1.0; // You can adjust this value

    let mut simulation_objects: Vec<SimulationObject> = Vec::new();

    simulation_objects.push(SimulationObject {
        shape: Circle { radius: 20.0 },
        selected: false,
        position: Vector2 {
            x: 700.0 as f32,
            y: 400.0 as f32,
        },
        body: RigidBody {
            mass: 200.0,
            velocity: Vector2 { x:  10.0, y: -20.0 },
            acceleration: Vector2 { x: 0.0, y: 0.0 },
        },
        material: Material {
            restitution: 1.0,
            colour: Color::BLUE,
        },
        history: History { trail: Vec::new() },
    });
    simulation_objects.push(SimulationObject {
        shape: Circle { radius: 33.0 },
        selected: false,
        position: Vector2 {
            x: 500.0 as f32,
            y: 400.0 as f32,
        },
        body: RigidBody {
            mass: 40.0,
            velocity: Vector2 { x: -50.0, y: 100.0 },
            acceleration: Vector2 { x: 0.0, y: 0.0 },
        },
        material: Material {
            restitution: 1.0,
            colour: Color::RED,
        },

        history: History { trail: Vec::new() },
    });

    // simulation_objects.push(SimulationObject {
    //     shape: Circle { radius: 10.0 },
    //     selected: false,
    //     position: Vector2 {
    //         x: 100.0 as f32,
    //         y: 200.0 as f32,
    //     },
    //     body: RigidBody {
    //         mass: 50.0,
    //         velocity: Vector2 { x: 0.0, y: 0.0 },
    //         acceleration: Vector2 { x: 0.0, y: 0.0 },
    //     },
    //     material: Material {
    //         restitution: 1.0,
    //         colour: Color::GREEN,
    //     },

    //     history: History { trail: Vec::new() },
    // });

    simulation_objects[selected_index].selected = true;

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();
    let timer_subsystem = sdl_context.timer().unwrap();
    let mut prev_mouse_pos = Vector2 { x: 0.0, y: 0.0 };

    let window = video_subsystem
        .window("SDL2 Demo", screen_width, screen_height)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();
    let mut event_pump = sdl_context.event_pump().unwrap();

    let mut last_ticks = timer_subsystem.ticks();

    let mut mouseDown = false;

    let mut predicted_positions = predict(&simulation_objects, &event_pump);
    let mut velocity_changed = false; // Add this flag
    let initial_simulation_objects = simulation_objects.clone();
    'running: loop {
        let current_ticks = timer_subsystem.ticks();
        let delta_time = (current_ticks - last_ticks) as f32 / 1000.0; // in seconds
        last_ticks = current_ticks;

        for event in event_pump.poll_iter() {
            match event {
                sdl2::event::Event::Quit { .. }
                | sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Escape),
                    ..
                } => {
                    break 'running;
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Tab),
                    ..
                } => {
                    simulation_objects[selected_index].selected = false;
                    selected_index = (selected_index + 1) % simulation_objects.len();
                    simulation_objects[selected_index].selected = true;
                }
                sdl2::event::Event::MouseButtonDown { x, y, .. } => {
                    mouseDown = true;
                    prev_mouse_pos = Vector2 {
                        x: x as f32,
                        y: y as f32,
                    };
                    for (i, object) in simulation_objects.iter().enumerate() {
                        if circles_intersect(
                            x as f32,
                            y as f32,
                            0.0,
                            object.position.x,
                            object.position.y,
                            object.shape.radius,
                        ) {
                            simulation_objects[selected_index].selected = false;
                            selected_index = i;
                            simulation_objects[selected_index].selected = true;
                            break;
                        }
                    }
                }

                sdl2::event::Event::MouseButtonUp { x, y, .. } => {
                    mouseDown = false;
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Space),
                    ..
                } => {
                    is_paused = !is_paused;
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::D),
                    ..
                } => {
                    single_step = true;
                }
                sdl2::event::Event::MouseMotion { x, y, .. } => {
                    if mouseDown {
                        let current_mouse_pos = Vector2 {
                            x: x as f32,
                            y: y as f32,
                        };
                        let dx = current_mouse_pos.x - prev_mouse_pos.x;
                        let dy = current_mouse_pos.y - prev_mouse_pos.y;

                        // Update all object positions
                        for object in simulation_objects.iter_mut() {
                            object.position.x += dx;
                            object.position.y += dy;
                            // Update all trail positions (assuming trails is a Vec<Vector2>)
                            for trail in &mut object.history.trail {
                                trail.x += dx;
                                trail.y += dy;
                            }
                        }

                        for p in predicted_positions.iter_mut() {
                            for q in p {
                                q.x += dx;
                                q.y += dy;
                            }
                        }

                        prev_mouse_pos = current_mouse_pos;
                    }
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Up),
                    ..
                } => {
                    simulation_objects[selected_index].body.velocity.y -= 10.0;
                    velocity_changed = true;
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Down),
                    ..
                } => {
                    simulation_objects[selected_index].body.velocity.y += 10.0;
                    velocity_changed = true;
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Left),
                    ..
                } => {
                    simulation_objects[selected_index].body.velocity.x -= 10.0;
                    velocity_changed = true;
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Right),
                    ..
                } => {
                    simulation_objects[selected_index].body.velocity.x += 10.0;
                    velocity_changed = true;
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::R),
                    ..
                } => {
                    simulation_objects = initial_simulation_objects.clone();
                    velocity_changed = true;
                }
                _ => {}
            }
        }

        if velocity_changed {
            predicted_positions = predict(&simulation_objects, &event_pump);
            velocity_changed = false; // Reset the flag
        }

        if !is_paused || single_step {
            single_step = false;

            physics::update_simulation(
                &mut simulation_objects,
                delta_time,
                selected_index,
                &event_pump,
            );
        }

        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();

        if is_paused {
            for object_index in 0..predicted_positions[0].len() {
                for i in 0..(predicted_positions.len() - 1) {
                    let start_pos = &predicted_positions[i][object_index];
                    let end_pos = &predicted_positions[i + 1][object_index];

                    canvas
                        .line(
                            start_pos.x as i16,
                            start_pos.y as i16,
                            end_pos.x as i16,
                            end_pos.y as i16,
                            Color::RGB(0, 255, 0), // Replace with your own color logic
                        )
                        .unwrap();
                }
            }
        }

        for object in simulation_objects.iter() {
            for point in &object.history.trail {
                canvas
                    .filled_circle(point.x as i16, point.y as i16, 2, object.material.colour)
                    .unwrap();
            }
        }

        for object in simulation_objects.iter_mut() {
            if object.selected {
                canvas
                    .filled_circle(
                        object.position.x as i16,
                        object.position.y as i16,
                        (object.shape.radius + 4.0) as i16,
                        Color::RGB(255, 255, 255),
                    )
                    .unwrap();
            }
            canvas
                .filled_circle(
                    object.position.x as i16,
                    object.position.y as i16,
                    object.shape.radius as i16,
                    object.material.colour,
                )
                .unwrap();
        }
        let texture_creator = canvas.texture_creator();
        let selected_object = &simulation_objects[selected_index];
        let info = format!(
            "A: ({:.2},{:.2}), V: ({:.2}, {:.2}), P: ({:.2}, {:.2})",
            selected_object.body.acceleration.x,
            selected_object.body.acceleration.y,
            selected_object.body.velocity.x,
            selected_object.body.velocity.y,
            selected_object.position.x,
            selected_object.position.y
        );
        render_text(&mut canvas, &texture_creator, &font, &info, 10, 10);
        canvas.present();
    }
}

fn check_collision_axis(
    current_index: usize,
    other_index: usize,
    future_positions: &[Vector2],
    objects: &[SimulationObject],
    axis: char,
) -> bool {
    let current_shape = objects[current_index].clone();
    let other_shape = objects[other_index].clone();

    let mut future_current_shape = current_shape;
    let mut future_other_shape = other_shape;

    match axis {
        'x' => {
            future_current_shape.position.x = future_positions[current_index].x;
            future_other_shape.position.x = future_positions[other_index].x;
        }
        'y' => {
            future_current_shape.position.y = future_positions[current_index].y;
            future_other_shape.position.y = future_positions[other_index].y;
        }
        _ => panic!("Invalid axis"),
    }

    have_intersection(future_current_shape, future_other_shape)
}

fn have_intersection(s_o1: SimulationObject, s_o2: SimulationObject) -> bool {
    circles_intersect(
        s_o1.position.x,
        s_o1.position.y,
        s_o1.shape.radius,
        s_o2.position.x,
        s_o2.position.y,
        s_o2.shape.radius,
    )
}

fn circles_intersect(x1: f32, y1: f32, r1: f32, x2: f32, y2: f32, r2: f32) -> bool {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let distance_squared = dx * dx + dy * dy;
    let radii_sum = r1 + r2;
    distance_squared <= radii_sum * radii_sum
}

fn apply_force(body: &mut RigidBody, force: Vector2) {
    body.acceleration.x += force.x / body.mass;
    body.acceleration.y += force.y / body.mass;
}

fn render_text(
    canvas: &mut Canvas<Window>,
    texture_creator: &sdl2::render::TextureCreator<sdl2::video::WindowContext>,
    font: &Font,
    text: &str,
    x: i32,
    y: i32,
) {
    let surface = font
        .render(text)
        .blended(Color::RGBA(255, 255, 255, 255))
        .unwrap();
    let texture = texture_creator
        .create_texture_from_surface(&surface)
        .unwrap();
    let TextureQuery { width, height, .. } = texture.query();
    let target = Rect::new(x, y, width, height);
    canvas.copy(&texture, None, Some(target)).unwrap();
}

fn predict(
    simulation_objects: &Vec<SimulationObject>,
    event_pump: &EventPump,
) -> Vec<Vec<Vector2>> {
    // Initialize an array to store the predicted positions
    let mut predicted_positions: Vec<Vec<Vector2>> = Vec::new();
    // Simulate the next N steps
    let steps = 2000;
    let mut temp_simulation_objects = simulation_objects.clone();
    for _ in 0..steps {
        let mut future_positions: Vec<Vector2> = Vec::new();

        // Create a deep copy of your simulation objects to use for prediction

        // Apply forces and update velocities and positions for each object
        physics::apply_forces(&mut temp_simulation_objects, 0, event_pump); // Assuming you have access to event_pump
        physics::update_velocities(&mut temp_simulation_objects, 0.1); // Assuming delta_time is known
        physics::handle_collisions(&mut temp_simulation_objects, 0.1);
        // Store the future positions
        for object in &temp_simulation_objects {
            future_positions.push(object.position);
        }

        // Store the future positions for this step
        predicted_positions.push(future_positions);
    }

    predicted_positions
}
