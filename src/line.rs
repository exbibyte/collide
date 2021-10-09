use crate::bound::Bound;
use crate::shape::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use num_traits::{Float, NumAssign};

use crate::bound_aabb::AxisAlignedBBox;
use lightmatrix::matrix::*;

#[derive(Debug, Clone)]
pub struct Line<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub _a: Matrix<T, 4, 1>,
    pub _b: Matrix<T, 4, 1>,
    pub _bound: AxisAlignedBBox<T>,
    pub _vicinity: T,
}

impl<T> Line<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn init(a: &[T], b: &[T]) -> Line<T> {
        assert!(a.len() == 3);
        assert!(b.len() == 3);

        let xs = vec![a[0], b[0]];
        let ys = vec![a[1], b[1]];
        let zs = vec![a[2], b[2]];

        use std::cmp::Ordering::*;

        let x_min = *xs
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
            .unwrap();
        let x_max = *xs
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
            .unwrap();

        let y_min = *ys
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
            .unwrap();
        let y_max = *ys
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
            .unwrap();

        let z_min = *zs
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
            .unwrap();
        let z_max = *zs
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
            .unwrap();

        Line {
            _a: Matrix::from([[a[0], a[1], a[2], T::one()]]).t(),
            _b: Matrix::from([[b[0], b[1], b[2], T::one()]]).t(),
            _bound: AxisAlignedBBox::new(
                ShapeType::Rect,
                &[x_min, y_min, z_min, x_max, y_max, z_max],
            ),
            _vicinity: T::epsilon(),
        }
    }
}

impl<T> AnyBase for Line<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for Line<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType {
        ShapeType::Line
    }
    fn get_bound(&self) -> &dyn Bound<T> {
        &self._bound
    }
    // this shall test for intersection of bounding shapes first before procedding to test intersection using algorithms of higher complexity
    fn get_intersect(&self, other: &dyn Shape<T>) -> (bool, Option<Matrix<T, 4, 1>>) {
        if !self.get_bound().intersect(other.get_bound()) {
            return (false, None);
        } else {
            match other.get_type() {
                // ShapeType::TriPrism => other.get_intersect(self),
                _ => {
                    unimplemented!();
                }
            }
        }
    }
    fn get_support(&self, _v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>> {
        unimplemented!();
    }
}

impl<T> Vicinity<T> for Line<T>
where
    T: NumAssign + Copy + Default + Float,
{
    fn set_vicinity(&mut self, epsilon: T) {
        self._vicinity = epsilon.abs();
    }
    fn within_vicinity(&self, a: T, b: T) -> bool {
        if a + self._vicinity >= b && a - self._vicinity <= b {
            true
        } else {
            false
        }
    }
}
