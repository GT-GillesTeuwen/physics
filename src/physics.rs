use std::collections::HashMap;
use std::f32::consts::PI;


use crate::quadtree;
use crate::quadtree::Quadtree;
use crate::quadtree::Rect;
// src/physics.rs
use crate::SimulationObject;
use crate::TargetPoint;
use crate::Vector2;
use crate::QUERY_RADIUS;
use rayon::prelude::*;

const G: f32 = 0.0;
pub fn apply_viscosity_forces(
    simulation_objects: &mut Vec<SimulationObject>,
    force_scale: &f32,
    smoothing_radius: &f32,
    delta_time: f32,
) {
    for i in 0..simulation_objects.len() {
        let viscosity_force =
            calculate_viscosity_force(&simulation_objects, smoothing_radius, i) * (*force_scale);
        if viscosity_force.magnitude() > 1000.0 {
            println!("[{}],{:?}", viscosity_force.magnitude(), viscosity_force);
        }
        simulation_objects[i].body.acceleration.x += viscosity_force.x * delta_time;
        simulation_objects[i].body.acceleration.y += viscosity_force.y * delta_time;
        simulation_objects[i].history.viscosity_force.push(Vector2{x:viscosity_force.x * delta_time,y:viscosity_force.y * delta_time})
    }
}

pub fn calculate_viscosity_force(
    simulation_objects: &Vec<SimulationObject>,
    smoothing_radius: &f32,
    particle_index: usize,
) -> Vector2 {
    let mut total_viscosity_force = Vector2 { x: 0.0, y: 0.0 };
    for i in 0..simulation_objects.len() {
        if i != particle_index {
            let dist = (simulation_objects.get(particle_index).unwrap().position
                - simulation_objects.get(i).unwrap().position)
                .magnitude();
            let influence = viscosity_smoothing(dist, *smoothing_radius);

            total_viscosity_force += (simulation_objects.get(i).unwrap().body.velocity
                - simulation_objects
                    .get(particle_index)
                    .unwrap()
                    .body
                    .velocity)
                * influence;
        }
    }

    return total_viscosity_force;
}

pub fn viscosity_smoothing(dist: f32, smoothing_radius: f32) -> f32 {
    if dist >= smoothing_radius {
        0.0
    } else {
        (smoothing_radius  - dist)
    }
}

pub fn apply_repulsion_forces(simulation_objects: &mut Vec<SimulationObject>, delta_time: f32) {
    let repulsion_coefficient = 100.0; // Coefficient for how strongly particles repel each other
    let minimum_distance = 2.0; // Minimum allowed distance before repulsion kicks in

    for i in 0..simulation_objects.len() {
        for j in i + 1..simulation_objects.len() {
            let obj1 = &simulation_objects[i];
            let obj2 = &simulation_objects[j];
            let dx = obj2.position.x - obj1.position.x;
            let dy = obj2.position.y - obj1.position.y;
            let distance_squared = dx * dx + dy * dy;
            let distance = distance_squared.sqrt();

            // Only apply repulsion if objects are closer than the minimum distance
            if distance < minimum_distance && distance > 0.0 {
                let overlap = minimum_distance - distance;
                let repulsion_force_magnitude = repulsion_coefficient * overlap / distance;

                let repulsion_force_x = repulsion_force_magnitude * dx;
                let repulsion_force_y = repulsion_force_magnitude * dy;

                // Apply repulsion force proportional to how much they overlap
                simulation_objects[i].body.acceleration.x -= repulsion_force_x * delta_time;
                simulation_objects[i].body.acceleration.y -= repulsion_force_y * delta_time;
                simulation_objects[j].body.acceleration.x += repulsion_force_x * delta_time;
                simulation_objects[j].body.acceleration.y += repulsion_force_y * delta_time;
            }
        }
    }
}

pub fn update_simulation_2<F>(
    objects: &mut Vec<SimulationObject>,
    pressure_force_scale: &f32,
    viscosity_force_scale: &f32,
    smoothing_radius: &f32,
    social_coefficient: &f32,
    cognitive_coefficient: &f32,
    momentum: bool,
    delta_time: f32,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    target_points: &Vec<TargetPoint>,
    optimisation_function: F,
) where
    F: Fn(f32, f32, &TargetPoint) -> f32,
{
    let query_radius = QUERY_RADIUS; // Define the radius within which to consider neighbors for the global best
    let neighbor_bests = calculate_neighbor_bests(objects, &(target_points[0]), optimisation_function, query_radius,(max_x,max_y));
    apply_neighbor_forces(objects, delta_time, neighbor_bests, social_coefficient, cognitive_coefficient);
    let gradients = calculate_all_gradients(objects, smoothing_radius);
    apply_pressure_forces(objects, &gradients, pressure_force_scale, delta_time);
    apply_viscosity_forces(objects, &viscosity_force_scale, smoothing_radius, delta_time);
    update_velocities(objects, delta_time, momentum);

    objects.iter_mut().for_each(|object| {
        update_object_position(object, delta_time);
    });
}


pub fn update_simulation<F>(
    objects: &mut Vec<SimulationObject>,
    pressure_force_scale: &f32,
    viscosity_force_scale: &f32,
    smoothing_radius: &f32,
    social_coefficient: &f32,
    cognitive_coefficient: &f32,
    momentum:bool,
    delta_time: f32,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    target_points: &Vec<TargetPoint>,
    optimisation_function: F,
)where
F: Fn(f32, f32, &TargetPoint) -> f32, {
    
    //let global_best = calculate_global_best(objects, target_points);
    let global_best = calculate_global_best_2(objects,&(target_points.get(0).unwrap()),optimisation_function);
    apply_forces(objects, delta_time, global_best,social_coefficient,cognitive_coefficient);
    //handle_wall_collisions(objects, min_x, max_x, min_y, max_y);
    let gradients = calculate_all_gradients(objects, smoothing_radius);
    //apply_repulsion_forces(objects, delta_time); // Prevent penetration
    apply_pressure_forces(objects, &gradients, pressure_force_scale, delta_time);
    apply_viscosity_forces(
        objects,
        &viscosity_force_scale,
        smoothing_radius,
        delta_time,
    );
    update_velocities(objects, delta_time,momentum);

    objects.iter_mut().for_each(|object| {
        update_object_position(object, delta_time);
    });
}

use rand::Rng;

pub fn apply_neighbor_forces(
    objects: &mut Vec<SimulationObject>,
    delta_time: f32,
    neighbor_bests: Vec<Vector2>,
    social_coefficient: &f32,
    cognitive_coefficient: &f32
) {
    objects.par_iter_mut().enumerate().for_each(|(index, obj)| {
        let mut rng = rand::thread_rng();
        // Vector to personal best
        let personal_vector = obj.best.position - obj.position;
        // Vector to neighbor best
        let neighbor_vector = neighbor_bests[index] - obj.position;

        // Random coefficients for cognitive and social components
        let rand_cognitive = rng.gen::<f32>();
        let rand_social = rng.gen::<f32>();

        // Apply cognitive and social components with random factors
        obj.body.acceleration += personal_vector * (*cognitive_coefficient * rand_cognitive) * delta_time;
        obj.body.acceleration += neighbor_vector * (*social_coefficient * rand_social) * delta_time;
        obj.history.best_force.push(personal_vector * (*cognitive_coefficient * rand_cognitive) * delta_time + neighbor_vector * (*social_coefficient * rand_social) * delta_time);
    });
}

pub fn apply_forces(objects: &mut Vec<SimulationObject>, delta_time: f32, global_best: Vector2, social_coefficient: &f32, cognitive_coefficient: &f32) {
  

    objects.par_iter_mut().for_each(|obj| {
        let mut rng = rand::thread_rng();
        // Vector to personal best
        let personal_vector = obj.best.position - obj.position;
        // Vector to global best
        let global_vector = global_best - obj.position;

        // println!("{:?}-{:?}={:?}",obj.best.position,obj.position,personal_vector);
        // println!("{:?}-{:?}={:?}\n",global_best,obj.position,global_vector);

        // Random coefficients for cognitive and social components
        let mut  rand_cognitive = rng.gen::<f32>();
        let mut rand_social = rng.gen::<f32>();

     
        // Apply cognitive and social components with random factors
        obj.body.acceleration += personal_vector * (*cognitive_coefficient * rand_cognitive) * delta_time;
        obj.body.acceleration += global_vector * (*social_coefficient * rand_social) * delta_time;
        obj.history.best_force.push(personal_vector * (*cognitive_coefficient * rand_cognitive) * delta_time + global_vector * (*social_coefficient * rand_social) * delta_time);
    });
}

pub fn calculate_global_best(
    objects: &mut Vec<SimulationObject>,
    target_points: &Vec<TargetPoint>,
) -> Vector2 {
    let mut global_best = objects[0].best.position; // Start with the first object's best position
    let mut min_distance = f32::MAX; // Initialize to the maximum float value

    for obj in objects.iter_mut() {
        let best_distance_for_object; // Initialize to the maximum float value for the current object

        // Determine the smallest distance to any target point, including penalties, for the current position
        let  current_best_distance = target_points
            .iter()
            .map(|target_point| {
                (obj.position - target_point.target).magnitude() + target_point.penalty
            })
            .fold(f32::MAX, f32::min);

        // Determine the smallest distance to any target point, including penalties, for the best position
        let best_distance = target_points
            .iter()
            .map(|target_point| {
                (obj.best.position - target_point.target).magnitude() + target_point.penalty
            })
            .fold(f32::MAX, f32::min);

        // If the current position's best calculated distance is better, update the object's best position
        if current_best_distance < best_distance {
            obj.best.position = obj.position;
            best_distance_for_object = current_best_distance;
        } else {
            best_distance_for_object = best_distance;
        }

        // If the best distance for this object is less than the global minimum distance, update the global best
        if best_distance_for_object < min_distance {
            min_distance = best_distance_for_object;
            global_best = obj.best.position;
        }
    }

    //println!("{:?} ({})", global_best, min_distance);
    global_best
}

pub fn calculate_global_best_2<F>(
    objects: &mut Vec<SimulationObject>,
    target_point: &TargetPoint,
    calculate_a_function: F,
) -> Vector2
where
    F: Fn(f32, f32, &TargetPoint) -> f32,
{
    let mut global_best = objects[0].best.position; // Start with the first object's best position
    let mut min_value = f32::MAX; // Initialize to the maximum float value

    for obj in objects.iter_mut() {
        let best_value_for_object; // Initialize to the maximum float value for the current object

        // Calculate the value of the function a(x, y) for the current position
        let current_value = calculate_a_function(obj.position.x, obj.position.y, target_point);

        // Calculate the value of the function a(x, y) for the best position
        let best_value = calculate_a_function(obj.best.position.x, obj.best.position.y, target_point);

        // If the current position's value is better, update the object's best position
        if current_value < best_value {
            obj.best.position = obj.position;
            best_value_for_object = current_value;
        } else {
            best_value_for_object = best_value;
        }

        // If the best value for this object is less than the global minimum value, update the global best
        if best_value_for_object < min_value {
            min_value = best_value_for_object;
            global_best = obj.best.position;
        }
    }
    //println!("{:?} ({})", global_best, min_value);
    global_best
}

pub fn calculate_neighbor_bests<F>(
    objects: &mut Vec<SimulationObject>,
    target_point: &TargetPoint,
    calculate_a_function: F,
    query_radius: f32,
    screen_dims: (f32,f32)
) -> Vec<Vector2>
where
    F: Fn(f32, f32, &TargetPoint) -> f32,
{
    let boundary = Rect { x: 0.0, y: 0.0, width: screen_dims.0, height: screen_dims.1 };
    let mut quadtree = Quadtree::new(boundary, 4);

    // Insert all objects into the quadtree
    for obj in objects.iter() {
        quadtree.insert(obj.clone());
    }

    let mut neighbor_bests = vec![Vector2 { x: 0.0, y: 0.0 }; objects.len()];

    for (index, obj) in objects.iter_mut().enumerate() {
        let query_range = Rect { x: obj.position.x, y: obj.position.y, width: query_radius, height: query_radius };
        let mut neighbors = Vec::new();
        quadtree.query(&query_range, &mut neighbors);

        let mut min_value = f32::MAX;
        let mut best_position = obj.position;

        for &neighbor in &neighbors {
            let best_value_for_neighbor;

            // Calculate the value of the function a(x, y) for the neighbor's position
            let current_value = calculate_a_function(neighbor.position.x, neighbor.position.y, target_point);

            // Calculate the value of the function a(x, y) for the best position
            let best_value = calculate_a_function(neighbor.best.position.x, neighbor.best.position.y, target_point);

            // If the neighbor's position's value is better, update the best position
            if current_value < best_value {
                best_value_for_neighbor = current_value;
                if best_value_for_neighbor < min_value {
                    min_value = best_value_for_neighbor;
                    best_position = neighbor.position;
                }
            } else {
                best_value_for_neighbor = best_value;
                if best_value_for_neighbor < min_value {
                    min_value = best_value_for_neighbor;
                    best_position = neighbor.best.position;
                }
            }
        }

        neighbor_bests[index] = best_position;
    }

    neighbor_bests
}





pub fn update_velocities(objects: &mut Vec<SimulationObject>, delta_time: f32,momentum:bool) {
    objects.par_iter_mut().for_each(|obj| {
        if momentum{
             // Update velocities based on acceleration
        obj.body.velocity.x += obj.body.acceleration.x * delta_time/1000.0;
        obj.body.velocity.y += obj.body.acceleration.y * delta_time/1000.0;
        }else{
             // Update velocities based on acceleration
        obj.body.velocity.x = obj.body.acceleration.x * delta_time/1.0;
        obj.body.velocity.y = obj.body.acceleration.y * delta_time/1.0;
        }
       

        //obj.body.velocity*=0.9995;
        // Reset acceleration after updating velocities
        obj.body.acceleration = Vector2 { x: 0.0, y: 0.0 };
    })
}

pub fn handle_wall_collisions(
    objects: &mut Vec<SimulationObject>,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
) {
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

pub fn calculate_total_energy(objects: &[SimulationObject]) -> f32 {
    let mut total_energy = 0.0;

    for obj in objects {
        // Calculate kinetic energy
        let speed_squared = obj.body.velocity.x.powf(2.0) + obj.body.velocity.y.powf(2.0);
        let kinetic_energy = 0.5 * obj.body.mass * speed_squared;

        //println!("{} {}",potential_energy,kinetic_energy);
        // Sum up the energies
        total_energy += kinetic_energy;
    }

    total_energy
}

fn update_object_position(obj: &mut SimulationObject, dt: f32) {
    //println!("Pre-update: Position ({}, {}) Velocity ({}, {})", obj.position.x, obj.position.y, obj.body.velocity.x, obj.body.velocity.y);
    obj.position.x += obj.body.velocity.x * dt;
    obj.position.y += obj.body.velocity.y * dt;

    //println!("Post-update: Position ({}, {}) Velocity ({}, {})", obj.position.x, obj.position.y, obj.body.velocity.x, obj.body.velocity.y);

    obj.history.trail.push(Vector2 {
        x: obj.position.x,
        y: obj.position.y,
    });
}

fn calculate_pressure_gradient(
    x: f32,
    y: f32,
    simulation_objects: &Vec<SimulationObject>,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    smoothing_radius: &f32,
) -> Vector2 {
    let delta = 1.0; // Small change in position for gradient calculation

    // Calculate density near the particle and near shifted positions
    let density_x_plus = calculate_density_with_walls(
        (x + delta),
        y,
        simulation_objects,
        min_x,
        max_x,
        min_y,
        max_y,
        *smoothing_radius,
    );
    let density_x_minus = calculate_density_with_walls(
        (x - delta),
        y,
        simulation_objects,
        min_x,
        max_x,
        min_y,
        max_y,
        *smoothing_radius,
    );
    let density_y_plus = calculate_density_with_walls(
        x,
        (y + delta),
        simulation_objects,
        min_x,
        max_x,
        min_y,
        max_y,
        *smoothing_radius,
    );
    let density_y_minus = calculate_density_with_walls(
        x,
        (y - delta),
        simulation_objects,
        min_x,
        max_x,
        min_y,
        max_y,
        *smoothing_radius,
    );
    // Compute gradients
    let gradient_x = (density_x_plus - density_x_minus) / (2.0 * delta);
    let gradient_y = (density_y_plus - density_y_minus) / (2.0 * delta);

    Vector2 {
        x: gradient_x,
        y: gradient_y,
    }
}

pub fn calculate_density_with_walls(
    x: f32,
    y: f32,
    simulation_objects: &Vec<SimulationObject>,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    smoothing_radius: f32,
) -> f32 {
    let mut total = calculate_density(x, y, simulation_objects, smoothing_radius);

    //Add additional density based on proximity to walls
    total += smoothing(max_x - x, 1.0); // Right wall
    total += smoothing(x - min_x, 1.0); // Left wall
    total += smoothing(max_y - y, 1.0); // Bottom wall
    total += smoothing(y - min_y, 1.0); // Top wall

    total
}

fn smoothing(dist: f32, smoothing_radius: f32) -> f32 {
    if dist >= smoothing_radius {
        0.0
    } else {
        // Increase 'density' as distance decreases
        let volume = PI * smoothing_radius.powi(4) / 6.0;
        (smoothing_radius - dist) * (smoothing_radius - dist) / volume
    }
}

fn calculate_all_gradients(
    simulation_objects: &Vec<SimulationObject>,
    smoothing_radius: &f32,
) -> Vec<Vector2> {
    simulation_objects
        .iter()
        .map(|obj| {
            calculate_pressure_gradient(
                obj.position.x,
                obj.position.y,
                simulation_objects,
                0.0,
                1920.0,
                0.0,
                1080.0,
                smoothing_radius,
            )
        })
        .collect()
}

fn apply_pressure_forces(
    simulation_objects: &mut Vec<SimulationObject>,
    gradients: &[Vector2],
    force_scale: &f32,
    delta_time: f32,
) {
    for (obj, gradient) in simulation_objects.iter_mut().zip(gradients.iter()) {
        // Calculate the magnitude of the gradient
       
        let gradient_magnitude = (gradient.x.powi(2) + gradient.y.powi(2)).sqrt();

        // Introduce a non-linear scaling factor for forces
        let force_magnitude = force_scale; // Using square of the magnitude

        // Normalize the gradient vector
        let normalized_gradient = if gradient_magnitude != 0.0 {
            Vector2 {
                x: gradient.x / gradient_magnitude,
                y: gradient.y / gradient_magnitude,
            }
        } else {
            Vector2 { x: 0.0, y: 0.0 }
        };

        // Apply scaled forces
        obj.body.acceleration.x += -normalized_gradient.x * force_magnitude * delta_time;
        obj.body.acceleration.y += -normalized_gradient.y * force_magnitude * delta_time;
        obj.history.pressure_force.push(Vector2 { x: -normalized_gradient.x * force_magnitude * delta_time, y: -normalized_gradient.y * force_magnitude * delta_time })
    }
}

fn calculate_density(
    x: f32,
    y: f32,
    simulation_objects: &Vec<SimulationObject>,
    smoothing_radius: f32,
) -> f32 {
    let mut total = 0.0;

    for obj in simulation_objects {
        if (x as f32- obj.position.x).abs()<0.00001 && (y as f32- obj.position.y).abs()<0.00001{
            continue;
        }
        let dx = x as f32 - obj.position.x;
        let dy = y as f32 - obj.position.y;
        let distance = (dx * dx + dy * dy).sqrt();
        total += smoothing(distance, smoothing_radius);
    }

    total
}
