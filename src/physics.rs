use sdl2::{keyboard::Scancode, EventPump};

use crate::{RigidBody, SimulationObject, Vector2};

const G: f32 = 10000.0;
const MAX_VELOCITY: f32 = 1000.0; // Choose an appropriate value

pub fn update_simulation(
    objects: &mut Vec<SimulationObject>,
    delta_time: f32,
    selected_index: usize,
    event_pump: &EventPump,
) {
    apply_forces(objects, selected_index, event_pump);
    update_velocities(objects, delta_time);
    handle_collisions(objects, delta_time);
    // ... other physics-related functions
}

pub fn apply_forces(
    simulation_objects: &mut Vec<SimulationObject>,
    selected_index: usize,
    event_pump: &EventPump,
) {
    let force_strength = 800.0;
    // Step 0: Apply forces to all objects to update accelerations
    let mut force = Vector2 { x: 0.0, y: 0.0 };

    let keys: sdl2::keyboard::KeyboardState = event_pump.keyboard_state();

    if keys.is_scancode_pressed(Scancode::W) {
        force.y = -force_strength * simulation_objects[selected_index].body.mass;
    }
    if keys.is_scancode_pressed(Scancode::A) {
        force.x = -force_strength * simulation_objects[selected_index].body.mass;
    }
    if keys.is_scancode_pressed(Scancode::S) {
        force.y = force_strength * simulation_objects[selected_index].body.mass;
    }
    if keys.is_scancode_pressed(Scancode::D) {
        force.x = force_strength * simulation_objects[selected_index].body.mass;
    }

    for object in simulation_objects.iter_mut() {
        // Reset acceleration
        object.body.acceleration = Vector2 { x: 0.0, y: 0.0 };
    }

    for i in 0..simulation_objects.len() {
        for j in i + 1..simulation_objects.len() {
            let dx = (simulation_objects[j].position.x - simulation_objects[i].position.x);
            let dy = (simulation_objects[j].position.y - simulation_objects[i].position.y);
            let distance_squared = dx * dx + dy * dy;
            let distance = distance_squared.sqrt();

            if distance == 0.0 {
                continue; // Skip if distance is zero to avoid division by zero
            }

            let force_magnitude =G
                * (simulation_objects[i].body.mass * simulation_objects[j].body.mass)
                / distance_squared;

            let force_x = force_magnitude * dx / distance;
            let force_y = force_magnitude * dy / distance;


            // Apply the force to both objects
            apply_force(
                &mut simulation_objects[i].body,
                Vector2 {
                    x: force_x,
                    y: force_y,
                },
            );
            apply_force(
                &mut simulation_objects[j].body,
                Vector2 {
                    x: -force_x,
                    y: -force_y,
                },
            );
        }
    }

    for object in simulation_objects.iter_mut() {
        //apply_force(&mut object.body, Vector2 { x: 0.0, y: 100.0 });
        if object.selected {
            apply_force(&mut object.body, force);
        }
    }
}

pub fn update_velocities(simulation_objects: &mut Vec<SimulationObject>, delta_time: f32) {
    // Step 1: Update velocities based on accelerations
    for object in simulation_objects.iter_mut() {
        object.body.velocity.x += object.body.acceleration.x * delta_time;
        object.body.velocity.y += object.body.acceleration.y * delta_time;
    }
}

pub fn handle_collisions(simulation_objects: &mut Vec<SimulationObject>, delta_time: f32) {
    let mut future_positions: Vec<Vector2> = Vec::new();
    for object in &mut *simulation_objects {
        let future_position = Vector2 {
            x: object.position.x + object.body.velocity.x * delta_time,
            y: object.position.y + object.body.velocity.y * delta_time,
        };
        future_positions.push(future_position);
    }

    // Step 3: Handle collisions and adjust velocities

    //// Step 3.1: Find all colliding pairs
    // Step 3: Handle collisions and adjust velocities
    let mut collisions: Vec<(usize, usize)> = Vec::new();
    for i in 0..simulation_objects.len() {
        for j in i + 1..simulation_objects.len() {
            if circles_intersect(
                future_positions[i].x,
                future_positions[i].y,
                simulation_objects[i].shape.radius,
                future_positions[j].x,
                future_positions[j].y,
                simulation_objects[j].shape.radius,
            ) {
                collisions.push((i, j));
            }
        }
    }

    // Handle collisions
    for (i, other_index) in collisions {
        // Handle collision in a general way, affecting both x and y components
        let dx = future_positions[i].x - future_positions[other_index].x;
        let dy = future_positions[i].y - future_positions[other_index].y;
        let distance = (dx * dx + dy * dy).sqrt();

        let nx = dx / distance;
        let ny = dy / distance;

        let dp_norm_1 =
            simulation_objects[i].body.velocity.x * nx + simulation_objects[i].body.velocity.y * ny;
        let dp_norm_2 = simulation_objects[other_index].body.velocity.x * nx
            + simulation_objects[other_index].body.velocity.y * ny;

        let m1 = simulation_objects[i].body.mass;
        let m2 = simulation_objects[other_index].body.mass;

        let m1_inv = 1.0 / m1;
        let m2_inv = 1.0 / m2;

        let optimized_p = ((2.0 * (dp_norm_1 - dp_norm_2)) / (m1_inv + m2_inv))
            * simulation_objects[i].material.restitution
            * simulation_objects[other_index].material.restitution;

        // Update velocities based on the collision
        simulation_objects[i].body.velocity.x -= optimized_p * m1_inv * nx;
        simulation_objects[i].body.velocity.y -= optimized_p * m1_inv * ny;

        simulation_objects[other_index].body.velocity.x += optimized_p * m2_inv * nx;
        simulation_objects[other_index].body.velocity.y += optimized_p * m2_inv * ny;

        // Update future_positions based on new velocities
        future_positions[i].x =
            simulation_objects[i].position.x + simulation_objects[i].body.velocity.x * delta_time;
        future_positions[i].y =
            simulation_objects[i].position.y + simulation_objects[i].body.velocity.y * delta_time;

        future_positions[other_index].x = simulation_objects[other_index].position.x
            + simulation_objects[other_index].body.velocity.x * delta_time;
        future_positions[other_index].y = simulation_objects[other_index].position.y
            + simulation_objects[other_index].body.velocity.y * delta_time;

        // Calculate separation vector
        let dx = future_positions[i].x - future_positions[other_index].x;
        let dy = future_positions[i].y - future_positions[other_index].y;
        let distance = (dx * dx + dy * dy).sqrt();
        let overlap = simulation_objects[i].shape.radius
            + simulation_objects[other_index].shape.radius
            - distance;

        if overlap > 0.0 {
            let separation_distance = overlap / 2.0;
            let angle = dy.atan2(dx);
            let separation_x = separation_distance * angle.cos();
            let separation_y = separation_distance * angle.sin();

            // Apply separation
            future_positions[i].x += separation_x;
            future_positions[i].y += separation_y;
            future_positions[other_index].x -= separation_x;
            future_positions[other_index].y -= separation_y;
        }
    }

    for (i, object) in simulation_objects.iter_mut().enumerate() {
        // // Handle X-axis wall collision
        // if future_positions[i].x - object.shape.radius < 0.0
        //     || future_positions[i].x + object.shape.radius > 1000 as f32
        // {
        //     object.body.velocity.x = -object.body.velocity.x;
        //     future_positions[i].x = object.position.x; // Reset to current position
        // }

        // // Handle Y-axis wall collision
        // if future_positions[i].y - object.shape.radius < 0.0
        //     || future_positions[i].y + object.shape.radius > 800 as f32
        // {
        //     object.body.velocity.y = -object.body.velocity.y;
        //     future_positions[i].y = object.position.y; // Reset to current position
        // }

        object.position = future_positions[i]; // Update to future position
    }

    update_pos_and_vel(simulation_objects, &mut future_positions)
}

pub fn update_pos_and_vel(
    simulation_objects: &mut Vec<SimulationObject>,
    future_positions: &mut Vec<Vector2>,
) {
    // Step 4: Update positions and velocities for all objects
    for (i, object) in simulation_objects.iter_mut().enumerate() {
        object.position = future_positions[i];

        if object.history.trail.len() > 10000 {
            // Keep only the last 100 positions
            object.history.trail.remove(0);
        }
        //object.history.trail.push(object.position);
    }

    for object in simulation_objects.iter_mut() {
        object.body.velocity.x = object.body.velocity.x.min(MAX_VELOCITY).max(-MAX_VELOCITY);
        object.body.velocity.y = object.body.velocity.y.min(MAX_VELOCITY).max(-MAX_VELOCITY);
        
    }
}

fn apply_force(body: &mut RigidBody, force: Vector2) {
    body.acceleration.x += force.x / body.mass;
    body.acceleration.y += force.y / body.mass;
}

fn circles_intersect(x1: f32, y1: f32, r1: f32, x2: f32, y2: f32, r2: f32) -> bool {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let distance_squared = dx * dx + dy * dy;
    let radii_sum = r1 + r2;
    distance_squared <= radii_sum * radii_sum
}
