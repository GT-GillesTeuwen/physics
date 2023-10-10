use sdl2::{keyboard::Scancode, EventPump};

use crate::{RigidBody, SimulationObject, Vector2};

const G: f32 = 10000.0;

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
            let dx = simulation_objects[j].position.x - simulation_objects[i].position.x;
            let dy = simulation_objects[j].position.y - simulation_objects[i].position.y;
            let distance_squared = dx * dx + dy * dy;
            let distance = distance_squared.sqrt();

            if distance == 0.0 {
                continue; // Skip if distance is zero to avoid division by zero
            }

            let force_magnitude = G
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
    let m = simulation_objects.len();
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
        let m1 = simulation_objects[i].body.mass;
        let m2 = simulation_objects[other_index].body.mass;

        // Handle collision in a general way, affecting both x and y components
        let dx = future_positions[i].x - future_positions[other_index].x;
        let dy = future_positions[i].y - future_positions[other_index].y;
        let distance = (dx * dx + dy * dy).sqrt();

        let nx = dx / distance;
        let ny = dy / distance;

        let tx = -ny;
        let ty = nx;

        let dpTan1 =
            simulation_objects[i].body.velocity.x * tx + simulation_objects[i].body.velocity.y * ty;
        let dpTan2 = simulation_objects[other_index].body.velocity.x * tx
            + simulation_objects[other_index].body.velocity.y * ty;

        let dpNorm1 =
            simulation_objects[i].body.velocity.x * nx + simulation_objects[i].body.velocity.y * ny;
        let dpNorm2 = simulation_objects[other_index].body.velocity.x * nx
            + simulation_objects[other_index].body.velocity.y * ny;

        let m1 = simulation_objects[i].body.mass;
        let m2 = simulation_objects[other_index].body.mass;

        let m1Inv = 1.0 / m1;
        let m2Inv = 1.0 / m2;

        let optimizedP = ((2.0 * (dpNorm1 - dpNorm2)) / (m1Inv + m2Inv))
            * simulation_objects[i].material.restitution
            * simulation_objects[other_index].material.restitution;

        // Update velocities based on the collision
        simulation_objects[i].body.velocity.x -= optimizedP * m1Inv * nx;
        simulation_objects[i].body.velocity.y -= optimizedP * m1Inv * ny;

        simulation_objects[other_index].body.velocity.x += optimizedP * m2Inv * nx;
        simulation_objects[other_index].body.velocity.y += optimizedP * m2Inv * ny;

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

    update_pos_and_vel(simulation_objects,&mut future_positions)

   
}

pub fn update_pos_and_vel(simulation_objects: &mut Vec<SimulationObject>,future_positions: &mut Vec<Vector2>){
     // Step 4: Update positions and velocities for all objects
     for (i, object) in simulation_objects.iter_mut().enumerate() {
        object.position = future_positions[i];

        if object.history.trail.len() > 10000 {
            // Keep only the last 100 positions
            object.history.trail.remove(0);
        }
        object.history.trail.push(object.position);
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

pub(crate) fn apply_forces_to_temp_object(
    temp_object: &mut SimulationObject,
    simulation_objects: &mut Vec<SimulationObject>,
    temp_simulation_objects: &mut Vec<SimulationObject>,
)->Vec<SimulationObject> {
    
    // Reset acceleration for all temporary objects
    for object in temp_simulation_objects.iter_mut() {
        object.body.acceleration = Vector2 { x: 0.0, y: 0.0 };
    }
    temp_object.body.acceleration = Vector2 { x: 0.0, y: 0.0 };

    // Calculate gravitational forces for all pairs
    for i in 0..temp_simulation_objects.len() {
        for j in 0..simulation_objects.len() {
            let dx = simulation_objects[j].position.x - temp_simulation_objects[i].position.x;
            let dy = simulation_objects[j].position.y - temp_simulation_objects[i].position.y;
            let distance_squared = dx * dx + dy * dy;
            let distance = distance_squared.sqrt();

            if distance == 0.0 {
                continue;
            }

            let force_magnitude = G * (temp_simulation_objects[i].body.mass * simulation_objects[j].body.mass) / distance_squared;
            let force_x = force_magnitude * dx / distance;
            let force_y = force_magnitude * dy / distance;

            // Apply the force to temp_simulation_objects[i]
            apply_force(
                &mut temp_simulation_objects[i].body,
                Vector2 {
                    x: force_x,
                    y: force_y,
                },
            );
        }
    }

    // Now apply forces to the temp_object based on all other objects
    for i in 0..temp_simulation_objects.len() {
        let dx = temp_simulation_objects[i].position.x - temp_object.position.x;
        let dy = temp_simulation_objects[i].position.y - temp_object.position.y;
        let distance_squared = dx * dx + dy * dy;
        let distance = distance_squared.sqrt();

        if distance == 0.0 {
            continue;
        }

        let force_magnitude = G * (temp_object.body.mass * temp_simulation_objects[i].body.mass) / distance_squared;
        let force_x = force_magnitude * dx / distance;
        let force_y = force_magnitude * dy / distance;

        // Apply the force to temp_object
        apply_force(
            &mut temp_object.body,
            Vector2 {
                x: force_x,
                y: force_y,
            },
        );
    }
    temp_simulation_objects.clone()
}

