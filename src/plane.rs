use crate::bound::Bound;
use crate::point::*;
use crate::shape::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use num_traits::{Float, NumAssign};

use crate::bound_aabb::AxisAlignedBBox;
use lightmatrix::matrix::*;

#[derive(Debug, Clone)]
pub struct Plane<T: NumAssign + Copy + Default + Float> {
    pub _offset: Matrix<T, 4, 1>,
    pub _normal: Matrix<T, 4, 1>,
    pub _bound: AxisAlignedBBox<T>,
    pub _vicinity: T,
}

impl<T> Plane<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn init(offset: &[T], normal: &[T]) -> Plane<T> {
        assert!(offset.len() == 3);
        assert!(normal.len() == 3);
        Plane {
            _offset: Matrix::from([[offset[0], offset[1], offset[2], T::one()]]).t(),
            _normal: Matrix::from([[normal[0], normal[1], normal[2], T::zero()]])
                .t()
                .normalize_l2(),
            _bound: AxisAlignedBBox::new(
                ShapeType::Plane,
                &[&offset[0..3], &normal[0..3]].concat(),
            ),
            _vicinity: T::epsilon(),
        }
    }
}

impl<T> AnyBase for Plane<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for Plane<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType {
        ShapeType::Plane
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
                ShapeType::Plane => {
                    unimplemented!();
                }
                ShapeType::Ray => {
                    //see Ray3 for ray plane intersection
                    return other.get_intersect(self);
                }
                ShapeType::Sphere => {
                    //see sphere for sphere plane intersection
                    return other.get_intersect(self);
                }
                ShapeType::Point => {
                    let other_point: &Point<T> = match other.as_any().downcast_ref::<Point<T>>() {
                        Some(b) => b,
                        None => {
                            panic!("cast to Point failed");
                        }
                    };
                    let b_off = other_point._ori;

                    let k = self._normal.inner(&self._offset);
                    let c = self._normal.inner(&b_off);
                    let d = k - c;
                    if !self.within_vicinity(d, T::zero()) {
                        return (false, None);
                    }
                    return (true, Some(b_off));
                }
                _ => {
                    unimplemented!();
                }
            }
        }
    }
    fn get_support(&self, _v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>> {
        None
    }
}

impl<T> Vicinity<T> for Plane<T>
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

#[test]
fn test_intersect_plane_point_0() {
    //plane point intersection
    {
        let a = Point::init(&[2f64, 1f64, 2f64]);
        let b = Plane::init(&[1f64, 1f64, 1f64], &[0f64, 1f64, 0f64]);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, Matrix::from([[2f64, 1f64, 2f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for plane point intersection"),
        }
    }
}
#[test]
fn test_intersect_plane_point_1() {
    //plane point intersection
    {
        let a = Point::init(&[2f64, 1f64, 2f64]);
        let b = Plane::init(&[1f64, 1f64, 1f64], &[0f64, 1f64, 0f64]);
        match b.get_intersect(&a) {
            (true, Some(loc)) => {
                assert_eq!(loc, Matrix::from([[2f64, 1f64, 2f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for plane point intersection"),
        }
    }
}
#[test]
fn test_intersect_plane_point_2() {
    //plane point no intersection
    {
        let a = Point::init(&[2f64, 1.05f64, 2f64]);
        let b = Plane::init(&[1f64, 1f64, 1f64], &[0f64, 1f64, 0f64]);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for plane point no intersection"),
        }
    }
}
#[test]
fn test_intersect_plane_point_3() {
    //plane point no intersection
    {
        let a = Point::init(&[2f64, 0.99f64, 2f64]);
        let b = Plane::init(&[1f64, 1f64, 1f64], &[0f64, 1f64, 0f64]);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for plane point no intersection"),
        }
    }
}
