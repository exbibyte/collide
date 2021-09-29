use crate::bound::Bound;
use crate::vicinity::Vicinity;
use core::any::Any;
use lightmatrix::matrix::*;
use num_traits::{Float, NumAssign};

pub trait AnyBase {
    fn as_any(&self) -> &dyn Any;
}

#[derive(Clone, Copy, Debug)]
pub enum ShapeType {
    //primitive shapes
    Point,
    Ray,
    Sphere,
    Plane,
    Trig,
    Box,
    Rect,
    TriPrism, //5 facets, 2 triangles, 3 rectangles
    Line,
    //todo
    Frustum,
    Complex, //custom shapes
}

pub trait Shape<T>: Vicinity<T> + AnyBase
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType;
    fn get_bound(&self) -> &dyn Bound<T>;
    //optionally returns a location of intersection of bounding shapes, preferrably closest of such locations
    fn get_intersect(&self, other: &dyn Shape<T>) -> (bool, Option<Matrix<T, 4, 1>>);
    //required for gjk intersection test
    fn get_support(&self, v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>>;
}
