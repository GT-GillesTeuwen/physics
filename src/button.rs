
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::{Point, Rect};
use sdl2::render::{Canvas, TextureCreator, TextureQuery};
use sdl2::ttf::Font;
use sdl2::video::{Window, WindowContext};
use sdl2::EventPump;

use crate::render_text;

pub struct Button {
    rect: Rect,
    label: String,
    on_click: Box<dyn FnMut()>,
}

impl Button {
    pub fn new(x: i32, y: i32, width: u32, height: u32, label: String, on_click: Box<dyn FnMut()>) -> Self {
        Button {
            rect: Rect::new(x, y, width, height),
            label,
            on_click,
        }
    }

    pub fn draw(&self, canvas: &mut Canvas<Window>, font: &Font, texture_creator: &TextureCreator<WindowContext>) {
        canvas.set_draw_color(Color::RGBA(200, 200, 200, 255));
        canvas.fill_rect(self.rect).expect("Could not fill rect for button");

        render_text(canvas, texture_creator, font, &self.label, self.rect.x() + 5, self.rect.y() + 5);
    }

    pub fn handle_event(&mut self, event: &sdl2::event::Event) {
        if let sdl2::event::Event::MouseButtonDown { x, y, mouse_btn: sdl2::mouse::MouseButton::Left, .. } = event {
            if self.rect.contains_point(Point::new(*x, *y)) {
                (self.on_click)();
            }
        }
    }
}
