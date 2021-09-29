use num_traits::{Float, NumAssign};

pub enum BoundType {
    AxisAlignBox,
    Sphere,
}

pub trait Bound<T>
where
    T: NumAssign + Copy + Default + Float,
{
    fn get_type(&self) -> BoundType;
    fn intersect(&self, other: &dyn Bound<T>) -> bool;
    fn get_shortest_separation(&self, other: &dyn Bound<T>) -> T;
    fn get_bound_data(&self) -> [T; 32];
    fn get_union(&mut self, bounds: &[&dyn Bound<T>]);
    fn get_centroid(&self) -> [T; 3];
}
