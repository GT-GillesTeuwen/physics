use events::handle_reset_to_initial_conditions;
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::pixels::Color;
use sdl2::rect::Rect;
use sdl2::render::{Canvas, TextureQuery};
use sdl2::ttf::Font;
use sdl2::video::Window;
use sdl2::EventPump;

mod events;
mod physics;
mod render;

enum EditableField {
    Mass,
    VelocityX,
    VelocityY,
    PositionX,
    PositionY,
    Radius,
    // Add more fields as needed
}

enum MenuState {
    Closed,
    Open(usize, EditableField, String), // Added EditableField here
}

enum AppState {
    Normal,
    Menu(MenuState),
}

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

    let screen_width = 1000;
    let screen_height = 800;
    let mut selected_index = 0;

    let mut simulation_objects: Vec<SimulationObject> = Vec::new();

    // for i in 0..500 {
    //     simulation_objects.push(SimulationObject {
    //         shape: Circle { radius: 7.0 },
    //         selected: false,
    //         position: Vector2 {
    //             x: (((i % 15) + 1) * 19) as f32,
    //             y: (((i / 15) + 1) * 19) as f32,
    //         },
    //         body: RigidBody {
    //             mass: 4.0,
    //             velocity: Vector2 { x: 0.0, y: 0.0 },
    //             acceleration: Vector2 { x: 0.0, y: 0.0 },
    //         },
    //         material: Material {
    //             restitution: 0.8,
    //             colour: Color::BLUE,
    //         },
    //         history: History { trail: Vec::new() },
    //     });
    // }
    simulation_objects.push(SimulationObject {
        shape: Circle { radius: 40.0 },
        selected: false,
        position: Vector2 {
            x: 500 as f32,
            y: 200 as f32,
        },
        body: RigidBody {
            mass: 100.0,
            velocity: Vector2 { x: 0.0, y: 0.0 },
            acceleration: Vector2 { x: 0.0, y: 0.0 },
        },
        material: Material {
            restitution: 0.8,
            colour: Color::RED,
        },
        history: History { trail: Vec::new() },
    });

    simulation_objects.push(SimulationObject {
        shape: Circle { radius: 40.0 },
        selected: false,
        position: Vector2 {
            x: 200 as f32,
            y: 700 as f32,
        },
        body: RigidBody {
            mass: 10.0,
            velocity: Vector2 { x: 0.0, y: 0.0 },
            acceleration: Vector2 { x: 0.0, y: 0.0 },
        },
        material: Material {
            restitution: 0.8,
            colour: Color::GREEN,
        },
        history: History { trail: Vec::new() },
    });
    simulation_objects[selected_index].selected = true;

    let mut camera_position = Vector2 { x: 0.0, y: 0.0 };

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

    let mut mouse_down = false;

    let mut app_state = AppState::Normal;
    let mut predicted_positions = predict(&simulation_objects, &event_pump);
    let mut repredict_path = false; // Add this flag
    let initial_simulation_objects = simulation_objects.clone();
    let mut just_toggled_menu = false;
    'running: loop {
        let selected_object = &simulation_objects[selected_index];
        let dx = screen_width as f32 / 2.0 - selected_object.position.x;
        let dy = screen_height as f32 / 2.0 - selected_object.position.y;
        let current_ticks = timer_subsystem.ticks();
        let mut delta_time = (current_ticks - last_ticks) as f32 / 1000.0; // in seconds
                                                                           //delta_time = 0.0075;
        last_ticks = current_ticks;

        for event in event_pump.poll_iter() {
            match &mut app_state {
                AppState::Normal => {
                    match event {
                        sdl2::event::Event::Quit { .. } => {
                            break 'running;
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::Tab),
                            ..
                        } => {
                            events::handle_tab_key(&mut simulation_objects, &mut selected_index);
                        }
                        sdl2::event::Event::MouseButtonDown { x, y, .. } => {
                            mouse_down = true;
                            events::handle_mouse_down(
                                x,
                                y,
                                &mut prev_mouse_pos,
                                &mut simulation_objects,
                                &mut selected_index,
                            );
                        }

                        sdl2::event::Event::MouseButtonUp { .. } => {
                            mouse_down = false;
                            events::handle_mouse_up();
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::Space),
                            ..
                        } => {
                            events::handle_space_key(&mut is_paused);
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::D),
                            ..
                        } => {
                            events::handle_single_step(&mut single_step);
                        }
                        sdl2::event::Event::MouseMotion { x, y, .. } => {
                            if mouse_down {
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
                            events::handle_velocity_change(
                                selected_index,
                                &mut simulation_objects,
                                0.0,
                                -5.0,
                            );
                            repredict_path = true;
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::Down),
                            ..
                        } => {
                            events::handle_velocity_change(
                                selected_index,
                                &mut simulation_objects,
                                0.0,
                                5.0,
                            );
                            repredict_path = true;
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::Left),
                            ..
                        } => {
                            events::handle_velocity_change(
                                selected_index,
                                &mut simulation_objects,
                                -5.0,
                                0.0,
                            );
                            repredict_path = true;
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::Right),
                            ..
                        } => {
                            events::handle_velocity_change(
                                selected_index,
                                &mut simulation_objects,
                                5.0,
                                0.0,
                            );
                            repredict_path = true;
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::R),
                            ..
                        } => {
                            handle_reset_to_initial_conditions(
                                &mut simulation_objects,
                                &initial_simulation_objects,
                            );
                            repredict_path = true;
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::M),
                            ..
                        } => {
                            just_toggled_menu = true;
                            match app_state {
                                AppState::Normal => {
                                    app_state = AppState::Menu(MenuState::Open(
                                        selected_index,
                                        EditableField::PositionX,
                                        simulation_objects[selected_index].position.x.to_string(),
                                    ));
                                }
                                AppState::Menu(_) => {
                                    app_state = AppState::Normal;
                                }
                            }
                        }
                        sdl2::event::Event::KeyDown {
                            keycode: Some(sdl2::keyboard::Keycode::C),
                            ..
                        } => {
                            let selected_object = &simulation_objects[selected_index];
                            camera_position.x =
                                (selected_object.position.x - screen_width as f32 / 2.0);
                            camera_position.y =
                                (selected_object.position.y - screen_height as f32 / 2.0);
                            repredict_path = true;
                        }
                        _ => {}
                    }
                    continue;
                }
                AppState::Menu(menu_state) => match menu_state {
                    MenuState::Closed => {}
                    MenuState::Open(index, field, input) => {
                        match event {
                            sdl2::event::Event::KeyDown {
                                keycode: Some(sdl2::keyboard::Keycode::Return),
                                ..
                            } => {
                                if let Ok(value) = input.parse::<f32>() {
                                    match field {
                                        EditableField::Mass => {
                                            simulation_objects[*index].body.mass = value;
                                        }
                                        EditableField::VelocityX => {
                                            simulation_objects[*index].body.velocity.x = value;
                                        }
                                        EditableField::VelocityY => todo!(),
                                        EditableField::PositionX => {
                                            simulation_objects[*index].position.x = value;
                                        }
                                        EditableField::PositionY => todo!(),
                                        EditableField::Radius => todo!(),
                                    }
                                }
                                app_state = AppState::Normal;
                                repredict_path = true;
                            }
                            sdl2::event::Event::KeyDown {
                                keycode: Some(sdl2::keyboard::Keycode::Backspace),
                                ..
                            } => {
                                if (input.len() > 0) {
                                    input.pop();
                                }
                            }
                            sdl2::event::Event::TextInput { text, .. } => {
                                if just_toggled_menu {
                                    just_toggled_menu = false;
                                    continue;
                                }
                                input.push_str(&text);
                            }
                            sdl2::event::Event::KeyDown {
                                keycode: Some(sdl2::keyboard::Keycode::Escape),
                                ..
                            } => {
                                app_state = AppState::Normal;
                            }
                            sdl2::event::Event::MouseButtonDown { x, y, .. } => {
                                mouse_down = true;
                                events::handle_mouse_down(
                                    x,
                                    y,
                                    &mut prev_mouse_pos,
                                    &mut simulation_objects,
                                    index, // Pass the selected index from the menu state
                                );
                            }
                            sdl2::event::Event::MouseButtonUp { .. } => {
                                mouse_down = false;
                                events::handle_mouse_up();
                            }
                            sdl2::event::Event::MouseMotion { x, y, .. } => {
                                if mouse_down {
                                    let current_mouse_pos = Vector2 {
                                        x: x as f32,
                                        y: y as f32,
                                    };
                                    let dx = current_mouse_pos.x - prev_mouse_pos.x;
                                    let dy = current_mouse_pos.y - prev_mouse_pos.y;

                                    // Update only the selected object's position
                                    simulation_objects[*index].position.x += dx;
                                    simulation_objects[*index].position.y += dy;

                                    prev_mouse_pos = current_mouse_pos;
                                    repredict_path = true;
                                }
                            }

                            _ => {}
                        }
                    }
                },
            }
        }

        if repredict_path {
            predicted_positions = predict(&simulation_objects, &event_pump);
            repredict_path = false; // Reset the flag
        }

        if !is_paused || single_step {
            single_step = false;

            for i in 1..=8{
            physics::update_simulation(
                &mut simulation_objects,
                (delta_time/8.0),
                selected_index,
                &event_pump,
            );
        }
            for object in simulation_objects.iter_mut() {
                object.position.x += dx;
                object.position.y += dy;
                // Update all trail positions (assuming trails is a Vec<Vector2>)
                for trail in &mut object.history.trail {
                    trail.x += dx;
                    trail.y += dy;
                }
            }

            for object_pos in &mut predicted_positions {
                for pos in object_pos{
                    pos.x+=dx;
                    pos.y+=dy;
                }
            }
        }

        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();

        if predicted_positions.len() > 0 && is_paused {
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
                            simulation_objects[object_index].material.colour, // Replace with your own color logic
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
            "A: ({:.2},{:.2}), V: ({:.2}, {:.2}), P: ({:.2}, {:.2}), M: {:.2}",
            selected_object.body.acceleration.x,
            selected_object.body.acceleration.y,
            selected_object.body.velocity.x,
            selected_object.body.velocity.y,
            selected_object.position.x,
            selected_object.position.y,
            selected_object.body.mass // Include mass here
        );
        render_text(&mut canvas, &texture_creator, &font, &info, 10, 10);

        match &app_state {
            AppState::Normal => {
                // ... (existing rendering for normal state)
            }
            AppState::Menu(menu_state) => {
                match menu_state {
                    MenuState::Closed => {
                        // Do nothing if the menu is closed
                    }
                    MenuState::Open(index, field, input) => {
                        // Draw your menu here based on the selected field
                        let menu_text = match field {
                            EditableField::Mass => format!("Edit Mass: {}", input),
                            EditableField::VelocityX => format!("Edit X-Velocity: {}", input),
                            EditableField::VelocityY => format!("Edit Y-Velocity: {}", input),
                            EditableField::PositionX => format!("Edit X-Position: {}", input),
                            EditableField::PositionY => format!("Edit Y-Position: {}", input),
                            EditableField::Radius => format!("Edit Radius: {}", input),
                        };
                        render_text(&mut canvas, &texture_creator, &font, &menu_text, 10, 50);
                    }
                }
            }
        }

        canvas.present();
    }
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
    let steps = 30000;
    let mut temp_simulation_objects = simulation_objects.clone();
    for _ in 0..steps {
        let mut future_positions: Vec<Vector2> = Vec::new();

        // Create a deep copy of your simulation objects to use for prediction

        // Apply forces and update velocities and positions for each object
        physics::apply_forces(&mut temp_simulation_objects, 0, event_pump); // Assuming you have access to event_pump
        physics::update_velocities(&mut temp_simulation_objects, 0.0075); // Assuming delta_time is known
        physics::handle_collisions(&mut temp_simulation_objects, 0.0075);
        // Store the future positions
        for object in &temp_simulation_objects {
            future_positions.push(object.position);
        }

        // Store the future positions for this step
        predicted_positions.push(future_positions);
    }

    predicted_positions
}
