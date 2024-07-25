use physics::calculate_total_energy;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use std::cell::RefCell;
use std::cmp;
use std::f32::consts::PI;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};
use std::rc::Rc; // Or use any other RNG from the `rand_pcg` crate

use button::Button;
use events::handle_reset_to_initial_conditions;
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;
use sdl2::render::{Canvas, TextureQuery};
use sdl2::ttf::Font;
use sdl2::video::Window;
use sdl2::EventPump;

mod button;
mod events;
mod physics;
mod quadtree;
mod render;

pub const QUERY_RADIUS:f32=150.0;

pub struct TargetPoint {
    target: Vector2,
    penalty: f32,
}

#[derive(Clone, Copy, Debug)]
pub struct Vector2 {
    x: f32,
    y: f32,
}
impl Vector2 {
    fn magnitude(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }
}

// Implement multiplication of Vector2 by f32
impl Mul<f32> for Vector2 {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Add<Vector2> for Vector2 {
    type Output = Self;

    fn add(self, rhs: Vector2) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl AddAssign<Vector2> for Vector2 {
    fn add_assign(&mut self, rhs: Vector2) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl Sub<Vector2> for Vector2 {
    type Output = Self;

    fn sub(self, rhs: Vector2) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl SubAssign<Vector2> for Vector2 {
    fn sub_assign(&mut self, rhs: Vector2) {
        self.x -= rhs.x;
        self.y -= rhs.y;
    }
}

// Implement multiplication assignment of Vector2 by f32
impl MulAssign<f32> for Vector2 {
    fn mul_assign(&mut self, rhs: f32) {
        self.x *= rhs;
        self.y *= rhs;
    }
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
    pressure_force: Vec<Vector2>,
    viscosity_force: Vec<Vector2>,
    best_force: Vec<Vector2>,
}
#[derive(Clone)]
pub struct Best {
    position: Vector2,
}

#[derive(Clone)]
pub struct SimulationObject {
    shape: Circle,
    selected: bool,
    position: Vector2,
    body: RigidBody,
    material: Material,
    history: History,
    best: Best,
}

impl SimulationObject {
    pub fn calculateEnergy(&self, max_y: f32, G: f32) -> f32 {
        (max_y - self.position.y - self.shape.radius) * G * &self.body.mass
            + 0.5
                * self.body.mass
                * (self.body.velocity.x.powf(2.0) + self.body.velocity.y.powf(2.0))
    }
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

    let font_path = "C:/Windows/Fonts/bahnschrift.ttf"; // Replace with your actual font path
    let font = ttf_context.load_font(font_path, 24).unwrap();

    let screen_width = 1920;
    let screen_height = 1080;

    let mut simulation_objects: Vec<SimulationObject> = Vec::new();
    let particle_radius = 2.0;
    let seed = [10, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; // A sample seed array
    let mut rng = Pcg64Mcg::from_seed(seed);

    for _ in 0..150 {
        // Assuming you want 200 particles as before (20 x 10)
        let x = rng.gen_range(0.0..screen_width as f32);
        let y = rng.gen_range(0.0..screen_height as f32);
        simulation_objects.push(SimulationObject {
            shape: Circle {
                radius: particle_radius,
            },
            selected: false,
            position: Vector2 { x: x, y: y },
            body: RigidBody {
                mass: 1.0,
                velocity: Vector2 { x: 0.0, y: 0.0 },
                acceleration: Vector2 { x: 0.0, y: 0.0 },
            },
            material: Material {
                restitution: 0.1,
                colour: Color {
                    r: 255,
                    g: 255,
                    b: 0,
                    a: 255,
                },
            },
            history: History {
                trail: Vec::new(),
                pressure_force: Vec::new(),
                viscosity_force: Vec::new(),
                best_force: Vec::new(),
            },
            best: Best {
                position: Vector2 { x: x, y: y },
            },
        });
    }

    let target_points = vec![
        TargetPoint {
            target: Vector2 {
                x: 960.0,
                y: 540.0,
            },
            penalty: 0.0,
        },
        // TargetPoint {
        //     target: Vector2 { x: 1200.0, y: 200.0 },
        //     penalty: 1.0,
        // },
        // TargetPoint {
        //     target: Vector2 { x: 350.0, y: 100.0 },
        //     penalty: 40.0,
        // },
        // TargetPoint {
        //     target: Vector2 { x: 800.0, y: 700.0 },
        //     penalty: 20.0,
        // },
        // TargetPoint {
        //     target: Vector2 { x: 1400.0, y: 600.0 },
        //     penalty: 20.0,
        // },
        // TargetPoint {
        //     target: Vector2 { x: 1600.0, y: 900.0 },
        //     penalty: 15.0,
        // },
    ];

    let pressure_force_scale = Rc::new(RefCell::new(500000.0));
    //let viscosity_force_scale = Rc::new(RefCell::new(0.010000));
    let viscosity_force_scale = Rc::new(RefCell::new(0.000000));
    let smoothing_radius = Rc::new(RefCell::new(100.0));
    let social_coefficient = Rc::new(RefCell::new(500.0));
    let cognitive_coefficient = Rc::new(RefCell::new(500.0));
    let momentum = Rc::new(RefCell::new(false));

    let momentum_toggle = momentum.clone();
    let mut toggle_momentum_button = Button::new(
        00,
        screen_height - 180,
        50,
        30,
        "!M".to_string(),
        Box::new(move || {
            let current_momentum = *momentum_toggle.borrow();
            *momentum_toggle.borrow_mut() = !current_momentum;
            println!("Toggled momentum: {}", momentum_toggle.borrow());
        }),
    );

    let social_coefficient_inc = social_coefficient.clone();
    let mut increase_social_coefficient_button = Button::new(
        60,
        screen_height - 150,
        50,
        30,
        "+SC".to_string(),
        Box::new(move || {
            *social_coefficient_inc.borrow_mut() *= 2.0;
            println!(
                "Increased social coefficient: {}",
                social_coefficient_inc.borrow()
            );
        }),
    );

    let social_coefficient_dec = social_coefficient.clone();
    let mut decrease_social_coefficient_button = Button::new(
        00,
        screen_height - 150,
        50,
        30,
        "-SC".to_string(),
        Box::new(move || {
            *social_coefficient_dec.borrow_mut() /= 2.0;
            println!(
                "Decreased social coefficient: {}",
                social_coefficient_dec.borrow()
            );
        }),
    );

    let cognitive_coefficient_inc = cognitive_coefficient.clone();
    let mut increase_cognitive_coefficient_button = Button::new(
        60,
        screen_height - 120,
        50,
        30,
        "+CC".to_string(),
        Box::new(move || {
            *cognitive_coefficient_inc.borrow_mut() *= 2.0;
            println!(
                "Increased cognitive coefficient: {}",
                cognitive_coefficient_inc.borrow()
            );
        }),
    );

    let cognitive_coefficient_dec = cognitive_coefficient.clone();
    let mut decrease_cognitive_coefficient_button = Button::new(
        0,
        screen_height - 120,
        50,
        30,
        "-CC".to_string(),
        Box::new(move || {
            *cognitive_coefficient_dec.borrow_mut() /= 2.0;
            println!(
                "Decreased cognitive coefficient: {}",
                cognitive_coefficient_dec.borrow()
            );
        }),
    );

    let smoothing_inc = smoothing_radius.clone();
    let mut increase_smoothing_radius_button = Button::new(
        60,
        screen_height - 90,
        50,
        30,
        "+S".to_string(),
        Box::new(move || {
            if *smoothing_inc.borrow() < 1.0 {
                *smoothing_inc.borrow_mut() += 0.1;
            } else {
                *smoothing_inc.borrow_mut() += 10.0;
            }

            println!("Increased smoothing radius: {}", smoothing_inc.borrow());
        }),
    );

    let smoothing_dec = smoothing_radius.clone();
    let mut decrease_smoothing_radius_button = Button::new(
        0,
        screen_height - 90,
        50,
        30,
        "-S".to_string(),
        Box::new(move || {
            if *smoothing_dec.borrow() <= 1.0 {
                *smoothing_dec.borrow_mut() -= 0.1;
            } else {
                *smoothing_dec.borrow_mut() -= 10.0;
            }
            println!("Decreased smoothing radius: {}", smoothing_dec.borrow());
        }),
    );

    let viscosity_inc = viscosity_force_scale.clone();
    let mut increase_viscosity_button = Button::new(
        60,
        screen_height - 60,
        50,
        30,
        "+V".to_string(),
        Box::new(move || {
            *viscosity_inc.borrow_mut() += 0.0000001;
            println!("Increased viscosity force: {}", viscosity_inc.borrow());
        }),
    );

    let viscosity_dec = viscosity_force_scale.clone();
    let mut decrease_viscosity_button = Button::new(
        0,
        screen_height - 60,
        50,
        30,
        "-V".to_string(),
        Box::new(move || {
            *viscosity_dec.borrow_mut() -= 0.0000001;
            println!("Decreased viscosity force: {}", viscosity_dec.borrow());
        }),
    );

    let pressure_inc = pressure_force_scale.clone();
    let mut increase_button = Button::new(
        60,
        screen_height - 30,
        50,
        30,
        "+P".to_string(),
        Box::new(move || {
            *pressure_inc.borrow_mut() += 1000.0;
            println!("Increased pressure force: {}", pressure_inc.borrow());
        }),
    );

    let pressure_dec = pressure_force_scale.clone();
    let mut decrease_button = Button::new(
        0,
        screen_height - 30,
        50,
        30,
        "-P".to_string(),
        Box::new(move || {
            if *pressure_dec.borrow() >= 10000.0 {
                *pressure_dec.borrow_mut() -= 10000.0;
            } else {
                *pressure_dec.borrow_mut() -= 1000.0;
            }
            println!("Decreased pressure force: {}", pressure_dec.borrow());
        }),
    );

    let pressure_zero = pressure_force_scale.clone();
    let mut zero_pressure_button = Button::new(
        120,
        screen_height - 30,
        50,
        30,
        "P=0".to_string(),
        Box::new(move || {
            *pressure_zero.borrow_mut() = 0.0;

            println!("0 pressure force: {}", pressure_zero.borrow());
        }),
    );

    let mut camera_position = Vector2 { x: 0.0, y: 0.0 };

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();
    let mut timer_subsystem = sdl_context.timer().unwrap();
    let mut prev_mouse_pos = Vector2 { x: 0.0, y: 0.0 };

    let window = video_subsystem
        .window("SDL2 Demo", screen_width, screen_height as u32)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();

    let texture_creator = canvas.texture_creator();
    let mut event_pump = sdl_context.event_pump().unwrap();

    let mut last_ticks = timer_subsystem.ticks();

    let fps = 60;
    let frame_delay = 1000 / fps;
    let updates_per_frame = 3.0;

    let mut ring_position: Option<(f32, f32)> = None;
    let mut push_radius: f32 = 100.0; // Assuming a default radius

    let mut is_dragging = false;
    let mut drag_position = (0.0, 0.0);

    let mut paused = true;
    let mut step_frame = false;

    let background_res = 15;
    let mut global_best:Vector2=Vector2 { x: 0.0, y: 0.0 };
    'running: loop {
        
        let start_ticks = timer_subsystem.ticks();

        // Clear the canvas
        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();

        for i in 0..=screen_width / background_res {
            for j in 0..=screen_height / background_res as i32 {
                let density = physics::calculate_density_with_walls(
                    (i * background_res) as f32,
                    (j as i32 * background_res as i32) as f32,
                    &simulation_objects,
                    0.0,
                    screen_width as f32,
                    0.0,
                    screen_height as f32,
                    *smoothing_radius.borrow(),
                );
                let color = density_to_colour(density, 0.2);
                canvas.set_draw_color(color);
                canvas
                    .fill_rect(sdl2::rect::Rect::new(
                        (i as i32 - 1) * background_res as i32,
                        (j as i32 - 1) * background_res as i32,
                        i * background_res,
                        (j * background_res as i32) as u32,
                    ))
                    .unwrap();
            }
        }

        for tp in &target_points {
            canvas
                .filled_circle(
                    tp.target.x as i16,
                    tp.target.y as i16,
                    4,
                    Color::RGB(
                        255 - (5.0 * tp.penalty) as u8,
                        255 - (5.0 * tp.penalty) as u8,
                        255 - (5.0 * tp.penalty) as u8,
                    ),
                )
                .unwrap();
        }

        increase_button.draw(&mut canvas, &font, &texture_creator);
        decrease_button.draw(&mut canvas, &font, &texture_creator);
        increase_viscosity_button.draw(&mut canvas, &font, &texture_creator);
        decrease_viscosity_button.draw(&mut canvas, &font, &texture_creator);
        increase_smoothing_radius_button.draw(&mut canvas, &font, &texture_creator);
        decrease_smoothing_radius_button.draw(&mut canvas, &font, &texture_creator);
        increase_social_coefficient_button.draw(&mut canvas, &font, &texture_creator);
        decrease_social_coefficient_button.draw(&mut canvas, &font, &texture_creator);
        increase_cognitive_coefficient_button.draw(&mut canvas, &font, &texture_creator);
        decrease_cognitive_coefficient_button.draw(&mut canvas, &font, &texture_creator);
        zero_pressure_button.draw(&mut canvas, &font, &texture_creator);
        toggle_momentum_button.draw(&mut canvas, &font, &texture_creator);
        // Event polling
        for event in event_pump.poll_iter() {
            match event {
                sdl2::event::Event::Quit { .. }
                | sdl2::event::Event::KeyDown {
                    keycode: Some(sdl2::keyboard::Keycode::Escape),
                    ..
                } => break 'running,
                sdl2::event::Event::MouseButtonDown {
                    mouse_btn: sdl2::mouse::MouseButton::Left,
                    x,
                    y,
                    ..
                } => {
                    increase_button.handle_event(&event);
                    decrease_button.handle_event(&event);
                    increase_viscosity_button.handle_event(&event);
                    decrease_viscosity_button.handle_event(&event);
                    increase_smoothing_radius_button.handle_event(&event);
                    decrease_smoothing_radius_button.handle_event(&event);
                    increase_social_coefficient_button.handle_event(&event);
                    decrease_social_coefficient_button.handle_event(&event);
                    increase_cognitive_coefficient_button.handle_event(&event);
                    decrease_cognitive_coefficient_button.handle_event(&event);
                    zero_pressure_button.handle_event(&event);
                    toggle_momentum_button.handle_event(&event);
                    let click_pos = Vector2 {
                        x: x as f32,
                        y: y as f32,
                    };
                    
                    for i in 0..simulation_objects.len() as i32 {
                        if (simulation_objects.get(i as usize).unwrap().position.x - click_pos.x)
                            .powi(2)
                            + (simulation_objects.get(i as usize).unwrap().position.y - click_pos.y)
                                .powi(2)
                            <= simulation_objects
                                .get(i as usize)
                                .unwrap()
                                .shape
                                .radius
                                .powi(2)
                        {
                            for o in &mut simulation_objects{
                                o.selected=false;
                            }
                            simulation_objects.get_mut(i as usize).unwrap().selected=true;
                            let obj=simulation_objects.get_mut(i as usize).unwrap();
                            println!(
                                "Forces on clicked particle at ({}, {}): {:?} {:?} {:?}",
                                obj.position.x,
                                obj.position.y,
                                obj
                                    .history
                                    .pressure_force
                                    .last(),
                                    obj
                                    .history
                                    .viscosity_force
                                    .last(),
                                    obj
                                    .history
                                    .best_force
                                    .last()
                            );
                        }
                    }
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(Keycode::Space),
                    ..
                } => {
                    paused = !paused; // Toggle pause
                }
                sdl2::event::Event::KeyDown {
                    keycode: Some(Keycode::N),
                    ..
                } => {
                    if paused {
                        step_frame = true; // Advance one frame
                    }
                }
                // Your existing mouse and motion event handling...
                _ => {
                    increase_button.handle_event(&event);
                    decrease_button.handle_event(&event);
                    increase_viscosity_button.handle_event(&event);
                    decrease_viscosity_button.handle_event(&event);
                    increase_smoothing_radius_button.handle_event(&event);
                    decrease_smoothing_radius_button.handle_event(&event);
                    increase_social_coefficient_button.handle_event(&event);
                    decrease_social_coefficient_button.handle_event(&event);
                    increase_cognitive_coefficient_button.handle_event(&event);
                    decrease_cognitive_coefficient_button.handle_event(&event);
                    zero_pressure_button.handle_event(&event);
                    toggle_momentum_button.handle_event(&event);
                }
            }
        }

        if !paused || step_frame {
            // Reset step_frame after advancing a frame
            if step_frame {
                step_frame = false;
            }


            if is_dragging {
                apply_force_to_objects(
                    &mut simulation_objects,
                    drag_position.0,
                    drag_position.1,
                    push_radius,
                    -100.0,
                );
            }
            // Update game logic and draw
            let dt = (timer_subsystem.ticks() - last_ticks) as f32 / 1000.0; // Convert milliseconds to seconds
            for i in 0..updates_per_frame as i32 {
                physics::update_simulation_2(
                    &mut simulation_objects,
                    &pressure_force_scale.borrow(),
                    &viscosity_force_scale.borrow(),
                    &smoothing_radius.borrow(),
                    &social_coefficient.borrow(),
                    &cognitive_coefficient.borrow(),
                    *momentum.borrow(),
                    dt / updates_per_frame,
                    0.0,
                    screen_width as f32,
                    0.0,
                    screen_height as f32,
                    &target_points,
                    ackley_function,
                );
            }
            if calculate_total_energy(&simulation_objects)<1.0 {
                break;
            }
            // Update game logic and draw
            // Update game logic and draw
            if is_dragging {
                canvas
                    .circle(
                        drag_position.0 as i16,
                        drag_position.1 as i16,
                        push_radius as i16,
                        Color::RGBA(255, 255, 255, 128),
                    )
                    .unwrap();
            }
        }

        let energy = physics::calculate_total_energy(&simulation_objects);
        let energy_text = format!("Total Energy: {:.7}", energy);
        let pressure_force_text = format!("Pressure Force: {:.2}", &pressure_force_scale.borrow());
        let viscosity_force_text =
            format!("Viscosity Force: {:.7}", &viscosity_force_scale.borrow());
        let smoothing_radius_text = format!("Smoothing Radius: {:.1}", &smoothing_radius.borrow());
        let cognitive_text = format!(
            "Cognitive Coefficient: {:.1}",
            &cognitive_coefficient.borrow()
        );
        let social_text = format!("Social Coefficient: {:.1}", &social_coefficient.borrow());
        let momentum_text = format!("Momentum: {}", &momentum.borrow());
        render_text(&mut canvas, &texture_creator, &font, &energy_text, 10, 10); // Position the text at (10, 10) on the canvas
        render_text(
            &mut canvas,
            &texture_creator,
            &font,
            &pressure_force_text,
            10,
            35,
        ); // Position the text at (10, 10) on the canvas
        render_text(
            &mut canvas,
            &texture_creator,
            &font,
            &viscosity_force_text,
            10,
            60,
        ); // Position the text at (10, 10) on the canvas
        render_text(
            &mut canvas,
            &texture_creator,
            &font,
            &smoothing_radius_text,
            10,
            85,
        ); // Position the text at (10, 10) on the canvas
        render_text(&mut canvas, &texture_creator, &font, &social_text, 10, 110); // Position the text at (10, 10) on the canvas
        render_text(
            &mut canvas,
            &texture_creator,
            &font,
            &cognitive_text,
            10,
            135,
        ); // Position the text at (10, 10) on the canvas
        render_text(
            &mut canvas,
            &texture_creator,
            &font,
            &momentum_text,
            10,
            160,
        ); // Position the text at (10, 10) on the canvas

        let neighbor_bests = physics::calculate_neighbor_bests(
            &mut simulation_objects,
            &(target_points[0]),
            ackley_function,
            QUERY_RADIUS,(screen_width as f32,screen_height as f32)
        );
        for n in neighbor_bests {
            canvas
                .filled_circle(
                    n.x as i16,
                    n.y as i16,
                    4.0 as i16,
                    Color::RGBA(100, 255, 100, 255),
                )
                .unwrap();
        }
        for (i,obj) in &mut simulation_objects.iter().enumerate() {
            draw_simulation_object(&mut canvas, obj);
            // let energy_text = format!("Particle Energy: {:.2}", obj.calculateEnergy(800.0, 9.8));
            // render_text(&mut canvas, &texture_creator, &font, &energy_text, 10, 30);
            // Position the text at (10, 10) on the canvas
            // canvas
            //     .circle(
            //         obj.position.x as i16,
            //         obj.position.y as i16,
            //         QUERY_RADIUS as i16,
            //         Color::RGBA(255, 100, 100, 75),
            //     )
            //     .unwrap();
            //     let a=smoothing_radius.borrow().clone() as i16;
            //     canvas
            //     .circle(
            //         obj.position.x as i16,
            //         obj.position.y as i16,
            //         a,
            //         Color::RGBA(51, 204, 255, 255),
            //     )
            //     .unwrap();
            if obj.selected{
                
                let neighbour_bests = physics::calculate_neighbor_bests(
                    &mut simulation_objects.clone(),
                    &(target_points[0]),
                    ackley_function,
                    QUERY_RADIUS,(screen_width as f32,screen_height as f32)
                );

                //println!("Particle Best: {:?}", neighbour_bests.get(i as usize));
                canvas
                    .filled_circle(
                        neighbour_bests.get(i as usize).unwrap().x as i16,
                        neighbour_bests.get(i as usize).unwrap().y as i16,
                        10.0 as i16,
                        Color::RGBA(100, 100, 255, 255),
                    )
                    .unwrap();
                //canvas.present();
            
            }
        }

        global_best = physics::calculate_global_best_2(
            &mut simulation_objects,
            &target_points.get(0).unwrap(),
            ackley_function,
        );
        canvas
            .filled_circle(
                global_best.x as i16,
                global_best.y as i16,
                4.0 as i16,
                Color::RGBA(255, 100, 100, 255),
            )
            .unwrap();

        // Update the canvas
        canvas.present();

        // Frame rate control and updating last_ticks for next loop
        let frame_time = timer_subsystem.ticks() - start_ticks;
        if frame_time < frame_delay {
            timer_subsystem.delay(frame_delay - frame_time);
        }
        last_ticks = start_ticks; // Update last_ticks for the next frame
    }

    println!("Optimal found: {:?}",global_best);
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

fn draw_simulation_object(canvas: &mut Canvas<Window>, simulation_object: &SimulationObject) {
    canvas.set_draw_color(simulation_object.material.colour);

    // Draw the circle representing the SimulationObject
    let center_x = simulation_object.position.x as i16;
    let center_y = simulation_object.position.y as i16;
    let radius = simulation_object.shape.radius as i16;

    // Using the filled_circle function from the gfx primitives to draw a filled circle
    canvas
        .filled_circle(
            center_x,
            center_y,
            radius,
            simulation_object.material.colour,
        )
        .unwrap();

    // Optionally, draw the trail if you need to visualize the object's historical positions
    for &pos in &simulation_object.history.trail {
        //canvas.filled_circle(pos.x as i16, pos.y as i16, 2, Color::RGBA(255, 255, 255, 128)).unwrap(); // Smaller and semi-transparent
    }
}

fn handle_mouse_click(
    objects: &mut [SimulationObject],
    mouse_x: f32,
    mouse_y: f32,
    push_radius: f32,
    push_strength: f32,
) -> (f32, f32) {
    for obj in objects.iter_mut() {
        let dx = obj.position.x - mouse_x;
        let dy = obj.position.y - mouse_y;
        let distance_squared = dx * dx + dy * dy;

        if distance_squared < push_radius * push_radius {
            let distance = distance_squared.sqrt();
            let push_factor = (push_radius - distance) / push_radius; // Normalize push factor
            let push_amount = push_strength * push_factor;

            let normalized_dx = dx / distance;
            let normalized_dy = dy / distance;

            // Apply force opposite to the direction of the cursor to object vector
            obj.body.velocity.x -= normalized_dx * push_amount;
            obj.body.velocity.y -= normalized_dy * push_amount;
        }
    }
    (mouse_x, mouse_y) // Return the mouse position
}

fn apply_force_to_objects(
    objects: &mut Vec<SimulationObject>,
    target_x: f32,
    target_y: f32,
    radius: f32,
    strength: f32,
) {
    for obj in objects.iter_mut() {
        let dx = target_x - obj.position.x;
        let dy = target_y - obj.position.y;
        let distance = (dx * dx + dy * dy).sqrt();

        if distance < radius {
            // Normalize direction
            let direction_x = dx / distance;
            let direction_y = dy / distance;

            // Apply a force that decreases with distance and considers object's velocity
            let force_magnitude = strength; // Avoid division by zero
            let damping = 0.2; // Damping factor to slow down as it approaches the target
            obj.body.velocity.x += (direction_x * force_magnitude - damping * obj.body.velocity.x);
            obj.body.velocity.y += (direction_y * force_magnitude - damping * obj.body.velocity.y);
        }
    }
}

fn density_to_colour(total: f32, target: f32) -> Color {
    let diff = (total * 1000.0 - target);
    // if total>0.0 {
    //     println!("{}-{}={}",total,target,diff);
    // }
    if diff > 0.0 {
        return Color::RGBA(cmp::min((diff * 100.0) as u8, 255) as u8, 0 as u8, 0, 255);
    } else if diff < 0.0 {
        return Color::RGBA(0 as u8, 0 as u8, cmp::min((diff * -100.0) as u8, 255), 255);
    }

    // Use the normalized influence to set the color intensity
    Color::RGBA(0, 255, 0, 200) // Keep blue constant for visibility
}

// Example function to calculate the value of a(x, y)
fn griewank_function(x: f32, y: f32, target_point: &TargetPoint) -> f32 {
    1.0 + (1.0 / 4000.0)
        * ((x - target_point.target.x).powi(2) + (y - target_point.target.y).powi(2))
        - ((x - target_point.target.x).cos() * ((y - target_point.target.y) / 2f32.sqrt()).cos())
}

// Define the 2D Rastrigin function
fn rastrigin_function(x: f32, y: f32, target_point: &TargetPoint) -> f32 {
    20.0 + (x - target_point.target.x).powi(2) + (y - target_point.target.y).powi(2)
        - 10.0
            * ((2.0 * std::f32::consts::PI * (x - target_point.target.x)).cos()
                + (2.0 * std::f32::consts::PI * (y - target_point.target.y)).cos())
}

fn deceptive_rastrigin(x: f32, y: f32, target_point: &TargetPoint) -> f32 {
    let a = 10.0;
    let b = 50.0; // Small noise factor to introduce irregularity

    let x_adj = (x - target_point.target.x) / 10.0;
    let y_adj = (y - target_point.target.y) / 10.0;

    // Standard Rastrigin computation
    let term1 = a * 2.0; // Scale factor for 2 dimensions
    let term2 = (x_adj * x_adj - a * (2.0 * PI * x_adj).cos())
        + (y_adj * y_adj - a * (2.0 * PI * y_adj).cos());

    // Adding noise to create deception
    let noise = b * ((x * y).sin() * 10.0).sin();

    term1 + term2 + noise
}

fn ackley_function(x: f32, y: f32, target_point: &TargetPoint) -> f32 {
    let a = 20.0;
    let b = 0.2;
    let c = 2.0 * std::f32::consts::PI;

    let x_adj = (x - target_point.target.x);
    let y_adj = (y - target_point.target.y);

    let term1 = -a * (-b * ((x_adj.powi(2) + y_adj.powi(2)) / 2.0).sqrt()).exp();
    let term2 = -(0.5 * ((c * x_adj).cos() + (c * y_adj).cos())).exp();

    term1 + term2 + a + std::f32::consts::E
}

fn custom_function_with_two_minima(x: f32, y: f32, target_point: &TargetPoint) -> f32 {
    let minima1 = (x - target_point.target.x+500.0).powi(2) + (y - target_point.target.y+300.0).powi(2)
        - 10.0 * ((2.0 * std::f32::consts::PI * (x - target_point.target.x+500.0)).cos()
        + (2.0 * std::f32::consts::PI * (y - target_point.target.y+300.0)).cos());
    
    let minima2 = (x - target_point.target.x-500.0).powi(2) + (y - target_point.target.y-300.0).powi(2)
        - 10.0 * ((2.0 * std::f32::consts::PI * (x - target_point.target.x-500.0)).cos()
        + (2.0 * std::f32::consts::PI * (y - target_point.target.y-300.0)).cos());

    // Return the minimum value of the two minima
    20.0 + minima1.min(minima2)
}

