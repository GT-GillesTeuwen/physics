use crate::SimulationObject;

// src/quadtree.rs
#[derive(Clone)]
pub struct Quadtree {
    boundary: Rect,
    capacity: usize,
    particles: Vec<SimulationObject>,
    divided: bool,
    northeast: Option<Box<Quadtree>>,
    northwest: Option<Box<Quadtree>>,
    southeast: Option<Box<Quadtree>>,
    southwest: Option<Box<Quadtree>>,
}

#[derive(Clone)]
pub struct Rect {
    pub x: f32,
    pub y: f32,
    pub width: f32,
    pub height: f32,
}

impl Rect {
    pub fn contains(&self, particle: &SimulationObject) -> bool {
        particle.position.x >= self.x - self.width &&
        particle.position.x < self.x + self.width &&
        particle.position.y >= self.y - self.height &&
        particle.position.y < self.y + self.height
    }

    pub fn intersects(&self, range: &Rect) -> bool {
        !(range.x - range.width > self.x + self.width ||
        range.x + range.width < self.x - self.width ||
        range.y - range.height > self.y + self.height ||
        range.y + range.height < self.y - self.height)
    }
}

impl Quadtree {
    pub fn new(boundary: Rect, capacity: usize) -> Self {
        Quadtree {
            boundary,
            capacity,
            particles: Vec::new(),
            divided: false,
            northeast: None,
            northwest: None,
            southeast: None,
            southwest: None,
        }
    }

    pub fn subdivide(&mut self) {
        let x = self.boundary.x;
        let y = self.boundary.y;
        let w = self.boundary.width / 2.0;
        let h = self.boundary.height / 2.0;

        let ne = Rect { x: x + w, y: y - h, width: w, height: h };
        let nw = Rect { x: x - w, y: y - h, width: w, height: h };
        let se = Rect { x: x + w, y: y + h, width: w, height: h };
        let sw = Rect { x: x - w, y: y + h, width: w, height: h };

        self.northeast = Some(Box::new(Quadtree::new(ne, self.capacity)));
        self.northwest = Some(Box::new(Quadtree::new(nw, self.capacity)));
        self.southeast = Some(Box::new(Quadtree::new(se, self.capacity)));
        self.southwest = Some(Box::new(Quadtree::new(sw, self.capacity)));

        self.divided = true;
    }

    pub fn insert(&mut self, particle: SimulationObject) -> bool {
        if !self.boundary.contains(&particle) {
            return false;
        }

        if self.particles.len() < self.capacity {
            self.particles.push(particle);
            return true;
        } else {
            if !self.divided {
                self.subdivide();
            }

            if self.northeast.as_mut().unwrap().insert(particle.clone()) {
                return true;
            } else if self.northwest.as_mut().unwrap().insert(particle.clone()) {
                return true;
            } else if self.southeast.as_mut().unwrap().insert(particle.clone()) {
                return true;
            } else if self.southwest.as_mut().unwrap().insert(particle) {
                return true;
            }
        }

        false
    }

    pub fn query<'a>(&'a self, range: &Rect, found: &mut Vec<&'a SimulationObject>) {
        if !self.boundary.intersects(range) {
            return;
        }

        for p in &self.particles {
            if range.contains(p) {
                found.push(p);
            }
        }

        if self.divided {
            self.northeast.as_ref().unwrap().query(range, found);
            self.northwest.as_ref().unwrap().query(range, found);
            self.southeast.as_ref().unwrap().query(range, found);
            self.southwest.as_ref().unwrap().query(range, found);
        }
    }
}
