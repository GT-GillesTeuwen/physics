use crate::Vector2;
use crate::SimulationObject;

pub fn handle_tab_key(simulation_objects: &mut Vec<SimulationObject>, selected_index: &mut usize) {
    simulation_objects[*selected_index].selected = false;
    *selected_index = (*selected_index + 1) % simulation_objects.len();
    simulation_objects[*selected_index].selected = true;
}

pub fn handle_mouse_down(x: i32, y: i32, prev_mouse_pos: &mut Vector2, simulation_objects: &mut Vec<SimulationObject>, selected_index: &mut usize) {
    *prev_mouse_pos = Vector2 { x: x as f32, y: y as f32 };
    for (i, object) in simulation_objects.iter().enumerate() {
        if circles_intersect(x as f32, y as f32, 0.0, object.position.x, object.position.y, object.shape.radius) {
            simulation_objects[*selected_index].selected = false;
            *selected_index = i;
            simulation_objects[*selected_index].selected = true;
            break;
        }
    }
}

pub fn handle_mouse_up() {
    // Your logic here, if any.
}

pub fn handle_space_key(is_paused: &mut bool) {
    *is_paused = !*is_paused;
}

pub fn handle_single_step(single_step: &mut bool) {
    *single_step = true;
}

pub fn handle_velocity_change(selected_index: usize, simulation_objects: &mut Vec<SimulationObject>, dx: f32, dy: f32) {
    simulation_objects[selected_index].body.velocity.x += dx;
    simulation_objects[selected_index].body.velocity.y += dy;
}

pub fn handle_reset_to_initial_conditions(simulation_objects: &mut Vec<SimulationObject>, initial_simulation_objects: &Vec<SimulationObject>) {
    *simulation_objects = initial_simulation_objects.clone();
}

// You can add more event logic functions as needed.

// Utility function for circle intersection
fn circles_intersect(x1: f32, y1: f32, r1: f32, x2: f32, y2: f32, r2: f32) -> bool {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let distance_squared = dx * dx + dy * dy;
    let radii_sum = r1 + r2;
    distance_squared <= radii_sum * radii_sum
}
