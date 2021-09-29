use num_traits::{Float, NumAssign};

pub trait Vicinity<T: NumAssign + Copy + Default + Float> {
    fn set_vicinity(&mut self, epsilon: T);
    fn within_vicinity(&self, a: T, b: T) -> bool;
}
