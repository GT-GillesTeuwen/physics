use std::collections::HashMap;

// src/physics.rs
use rayon::prelude::*;
use crate::SimulationObject;
use crate::Vector2;

const G:f32=0.0;







pub fn update_simulation(
    objects: &mut Vec<SimulationObject>,
    delta_time: f32,
    min_x: f32, max_x: f32, min_y: f32, max_y: f32
) {
   
    apply_forces(objects);  
    handle_wall_collisions(objects, min_x, max_x, min_y, max_y);
    let gradients = calculate_all_gradients(objects);
    apply_pressure_forces(objects, &gradients);
    update_velocities(objects, delta_time);

    objects.iter_mut().for_each(|object|{
        update_object_position(object, delta_time);
    });
}

pub fn apply_forces(objects: &mut Vec<SimulationObject>) {
    // Apply gravitational forces, wind, etc
    objects.par_iter_mut().for_each(|obj|
    {
        obj.body.acceleration.y += G; // Gravity
    })
}

pub fn update_velocities(objects: &mut Vec<SimulationObject>, delta_time: f32) {
    objects.par_iter_mut().for_each(|obj|
        {
             // Update velocities based on acceleration
        obj.body.velocity.x += obj.body.acceleration.x * delta_time;
        obj.body.velocity.y += obj.body.acceleration.y * delta_time;

        // Reset acceleration after updating velocities
        obj.body.acceleration = Vector2 { x: 0.0, y: 0.0 };
        })
    
}

pub fn handle_collisions(objects: &mut [SimulationObject], delta_time: f32) {
    let len = objects.len();
    for i in 0..len {
        let (left, right) = objects.split_at_mut(i + 1);
        if let Some(obj1) = left.last_mut() {
            for obj2 in right {
                if let Some(collision) = check_collision(obj1, obj2) {
                    resolve_collision(obj1, obj2, collision, delta_time);
                }
            }
        }
    }
}

fn check_collision(obj1: &SimulationObject, obj2: &SimulationObject) -> Option<(Vector2, f32)> {
    let dx = obj2.position.x - obj1.position.x;
    let dy = obj2.position.y - obj1.position.y;
    let distance_squared = dx * dx + dy * dy;

    let radius_sum = obj1.shape.radius + obj2.shape.radius;

    // Check if the square of the sum of the radii is greater than the square of the distance
    if distance_squared <= radius_sum * radius_sum {
        let distance = distance_squared.sqrt();

        if distance != 0.0 {
            // Normal vector of the collision (normalized vector between centers)
            let normal = Vector2 {
                x: dx / distance,
                y: dy / distance,
            };
            // Penetration depth (how much the circles overlap)
            let penetration = radius_sum - distance;
            Some((normal, penetration))
        } else {
            // Special case: circles are exactly on top of each other
            // Choose an arbitrary direction if centers are exactly the same
            Some((Vector2 { x: 1.0, y: 0.0 }, radius_sum))
        }
    } else {
        None
    }
}

fn resolve_collision(obj1: &mut SimulationObject, obj2: &mut SimulationObject, collision: (Vector2, f32), _delta_time: f32) {
    let (normal, penetration) = collision;
    let restitution = (obj1.material.restitution + obj2.material.restitution) / 2.0;

    // Calculate relative velocity
    let rel_velocity = Vector2 {
        x: obj2.body.velocity.x - obj1.body.velocity.x,
        y: obj2.body.velocity.y - obj1.body.velocity.y,
    };

    // Calculate velocity along the normal direction
    let vel_along_normal = rel_velocity.x * normal.x + rel_velocity.y * normal.y;

    // Do not resolve if velocities are separating
    if vel_along_normal > 0.0 {
        return;
    }

    // Calculate impulse scalar
    let impulse = -(1.0 + restitution) * vel_along_normal / (1.0 / obj1.body.mass + 1.0 / obj2.body.mass);

    // Apply impulse to the velocity components of both objects
    let impulse_x = impulse * normal.x;
    let impulse_y = impulse * normal.y;

    obj1.body.velocity.x -= impulse_x / obj1.body.mass;
    obj1.body.velocity.y -= impulse_y / obj1.body.mass;
    obj2.body.velocity.x += impulse_x / obj2.body.mass;
    obj2.body.velocity.y += impulse_y / obj2.body.mass;

    // Positional correction to remove penetration
    const POSITIONAL_CORRECTION: f32 = 0.5; // Use half for each object to push them out of penetration
    let correction = normal * penetration * POSITIONAL_CORRECTION;
    obj1.position.x -= correction.x;
    obj1.position.y -= correction.y;
    obj2.position.x += correction.x;
    obj2.position.y += correction.y;
}



pub fn handle_wall_collisions(objects: &mut Vec<SimulationObject>, min_x: f32, max_x: f32, min_y: f32, max_y: f32) {
    for obj in objects.iter_mut() {
        let radius = obj.shape.radius;
        // Horizontal walls
        if obj.position.x - radius <= min_x {
            obj.position.x = min_x + radius;
            obj.body.velocity.x = -obj.body.velocity.x * obj.material.restitution;
        } else if obj.position.x + radius >= max_x {
            obj.position.x = max_x - radius;
            obj.body.velocity.x = -obj.body.velocity.x * obj.material.restitution;
        }

        // Vertical walls
        if obj.position.y - radius <= min_y {
            obj.position.y = min_y + radius;
            obj.body.velocity.y = -obj.body.velocity.y * obj.material.restitution;
        } else if obj.position.y + radius >= max_y {
            //println!("{} {}",obj.position.y,obj.body.velocity.y);
            obj.position.y = max_y - radius;
            obj.body.velocity.y = -obj.body.velocity.y * obj.material.restitution;
            //println!("{} {}",obj.position.y,obj.body.velocity.y);
        }
    }
}

pub fn calculate_total_energy(objects: &[SimulationObject], max_y: f32) -> f32 {
    let mut total_energy = 0.0;

    for obj in objects {
        // Calculate kinetic energy
        let speed_squared = obj.body.velocity.x.powf(2.0) + obj.body.velocity.y.powf(2.0);
        let kinetic_energy = 0.5 * obj.body.mass * speed_squared;

        // Calculate potential energy, adjusted for graphical coordinate system where y increases downwards
        let height = max_y - obj.position.y-obj.shape.radius;
        let potential_energy = obj.body.mass * G * height;

        //println!("{} {}",potential_energy,kinetic_energy);
        // Sum up the energies
        total_energy += kinetic_energy + potential_energy;
    }

    total_energy
}




fn update_object_position(obj: &mut SimulationObject, dt: f32) {
    //println!("Pre-update: Position ({}, {}) Velocity ({}, {})", obj.position.x, obj.position.y, obj.body.velocity.x, obj.body.velocity.y);
    obj.position.x += obj.body.velocity.x * dt;
    obj.position.y += obj.body.velocity.y * dt;
    
    //println!("Post-update: Position ({}, {}) Velocity ({}, {})", obj.position.x, obj.position.y, obj.body.velocity.x, obj.body.velocity.y);

    obj.history.trail.push(Vector2 { x: obj.position.x, y: obj.position.y });
}

fn check_future_collision(obj1: &SimulationObject, obj2: &SimulationObject, dt: f32) -> Option<(Vector2, f32)> {
    // Predict future positions
    let future_position1 = Vector2 {
        x: obj1.position.x + obj1.body.velocity.x * dt,
        y: obj1.position.y + obj1.body.velocity.y * dt,
    };
    let future_position2 = Vector2 {
        x: obj2.position.x + obj2.body.velocity.x * dt,
        y: obj2.position.y + obj2.body.velocity.y * dt,
    };

    let dx = future_position2.x - future_position1.x;
    let dy = future_position2.y - future_position1.y;
    let distance_squared = dx * dx + dy * dy;
    let radius_sum = obj1.shape.radius + obj2.shape.radius;

    // Check if the square of the sum of the radii is greater than the square of the distance
    if distance_squared <= radius_sum * radius_sum {
        let distance = distance_squared.sqrt();

        // Normal vector of the collision (normalized vector between centers)
        let normal = Vector2 {
            x: dx / distance,
            y: dy / distance,
        };
        // Penetration depth (how much the circles overlap)
        let penetration = radius_sum - distance;
        Some((normal, penetration))
    } else {
        None
    }
}

fn handle_future_collisions(objects: &mut [SimulationObject], delta_time: f32) {
    let len = objects.len();
    for i in 0..len {
        let (left, right) = objects.split_at_mut(i + 1);
        if let Some(obj1) = left.last_mut() {
            for obj2 in right {
                if let Some(collision) = check_future_collision(obj1, obj2, delta_time) {
                    resolve_collision(obj1, obj2, collision, delta_time);
                }
            }
        }
    }
}


fn calculate_pressure_gradient(x: f32, y: f32, simulation_objects: &Vec<SimulationObject>, min_x: f32, max_x: f32, min_y: f32, max_y: f32) -> Vector2 {
    let delta = 50.0; // Small change in position for gradient calculation
    let smoothing_radius = 1000.0; // Example smoothing radius for density calculation

    // Calculate density near the particle and near shifted positions
    let density_x_plus = calculate_density_with_walls((x + delta), y, simulation_objects, min_x, max_x, min_y, max_y, smoothing_radius);
    let density_x_minus = calculate_density_with_walls((x - delta), y, simulation_objects, min_x, max_x, min_y, max_y, smoothing_radius);
    let density_y_plus = calculate_density_with_walls(x, (y + delta), simulation_objects, min_x, max_x, min_y, max_y, smoothing_radius);
    let density_y_minus = calculate_density_with_walls(x, (y - delta), simulation_objects, min_x, max_x, min_y, max_y, smoothing_radius);

    // Compute gradients
    let gradient_x = (density_x_plus - density_x_minus) / (2.0 * delta);
    let gradient_y = (density_y_plus - density_y_minus) / (2.0 * delta);

    Vector2 { x: gradient_x, y: gradient_y }
}

fn calculate_density_with_walls(x: f32, y: f32, simulation_objects: &Vec<SimulationObject>, min_x: f32, max_x: f32, min_y: f32, max_y: f32, smoothing_radius: f32) -> f32 {
    let mut total = 10.0*calculate_density(x as u32, y as u32, simulation_objects);

    // Add additional density based on proximity to walls
    total += 20.0*smoothing(max_x - x, smoothing_radius); // Right wall
    total += 20.0*smoothing(x - min_x, smoothing_radius); // Left wall
    total += 20.0*smoothing(max_y - y, smoothing_radius); // Bottom wall
    total += 20.0*smoothing(y - min_y, smoothing_radius); // Top wall

    total
}

fn smoothing(dist: f32, smoothing_radius: f32) -> f32 {
    if dist >= smoothing_radius {
        0.0
    } else {
        // Increase 'density' as distance decreases
        5.0*(smoothing_radius - dist)
    }
}


fn calculate_all_gradients(simulation_objects: &Vec<SimulationObject>) -> Vec<Vector2> {
    simulation_objects.iter().map(|obj| {
        calculate_pressure_gradient(obj.position.x, obj.position.y, simulation_objects,0.0,1000.0,0.0,800.0)
    }).collect()
}


fn apply_pressure_forces(simulation_objects: &mut Vec<SimulationObject>, gradients: &[Vector2]) {
    for (obj, gradient) in simulation_objects.iter_mut().zip(gradients.iter()) {
        // Assuming negative gradient direction for force
        obj.body.velocity.x = -gradient.x;
        obj.body.velocity.y = -gradient.y;
    }
}

fn calculate_density(x: u32, y: u32, simulation_objects: &Vec<SimulationObject>) -> f32 {
    
    let mut total=0.0;

    for obj in simulation_objects{
        let dx = x as f32 - obj.position.x-obj.shape.radius;
        let dy = y as f32 - obj.position.y-obj.shape.radius;
        let distance = (dx * dx + dy * dy).sqrt();
        total+=smoothing(distance, 100.0);
    }

   total
}

