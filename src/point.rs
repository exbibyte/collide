use crate::bound::Bound;
use crate::shape::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use num_traits::{Float, NumAssign};

use crate::bound_aabb::AxisAlignedBBox;
use lightmatrix::matrix::*;

#[derive(Debug, Clone)]
pub struct Point<T: NumAssign + Copy + Default + Float> {
    pub _ori: Matrix<T, 4, 1>,
    pub _bound: AxisAlignedBBox<T>,
    pub _vicinity: T,
}

impl<T> Point<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn init(origin: &[T]) -> Point<T> {
        assert!(origin.len() == 3);
        Point {
            _ori: Matrix::from([[origin[0], origin[1], origin[2], T::one()]]).t(),
            _bound: AxisAlignedBBox::new(ShapeType::Point, &origin[0..3]),
            _vicinity: T::from(1e-7).unwrap(),
        }
    }
}

impl<T> AnyBase for Point<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for Point<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType {
        ShapeType::Point
    }
    fn get_bound(&self) -> &dyn Bound<T> {
        &self._bound
    }
    // this shall test for intersection of bounding shapes first before procedding to test intersection using algorithms of higher complexity
    fn get_intersect(&self, other: &dyn Shape<T>) -> (bool, Option<Matrix<T, 4, 1>>) {
        if !self.get_bound().intersect(other.get_bound()) {
            (false, None)
        } else {
            match other.get_type() {
                ShapeType::Point => {
                    let other_point: &Point<T> = match other.as_any().downcast_ref::<Self>() {
                        Some(b) => b,
                        None => {
                            panic!("cast to Point failed");
                        }
                    };
                    let mut test = true;
                    for i in 0..3 {
                        test &= self.within_vicinity(self._ori[[i, 0]], other_point._ori[[i, 0]]);
                    }
                    if test {
                        (true, Some(self._ori))
                    } else {
                        (false, None)
                    }
                }
                ShapeType::Ray => {
                    //see Ray for ray point intersection
                    other.get_intersect(self)
                }
                ShapeType::Sphere => {
                    //see sphere for sphere point intersection
                    other.get_intersect(self)
                }
                // ShapeType::Plane => {
                //     //see plane for plane point intersection
                //     other.get_intersect(self)
                // }
                // ShapeType::Box => {
                //     //see recbox for box point intersection
                //     other.get_intersect(self)
                // }
                // ShapeType::TriPrism => {
                //     //see tri prism for intersection
                //     other.get_intersect(self)
                // }
                _ => {
                    unimplemented!();
                }
            }
        }
    }
    fn get_support(&self, _v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>> {
        Some(self._ori)
    }
}

impl<T> Vicinity<T> for Point<T>
where
    T: NumAssign + Copy + Default + Float,
{
    fn set_vicinity(&mut self, epsilon: T) {
        self._vicinity = epsilon.abs();
    }
    fn within_vicinity(&self, a: T, b: T) -> bool {
        a + self._vicinity >= b && a - self._vicinity <= b
    }
}

#[test]
fn test_intersect_point_point() {
    //point point intersection
    {
        let a = Point::init(&[25f64, 5f64, 5f64]);
        let b = Point::init(&[25f64, 5f64, 5f64]);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == b._ori);
            }
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
    //point point no intersection
    {
        let a = Point::init(&[25f64, 5f64, 5f64]);
        let b = Point::init(&[25.1f64, 5f64, 5f64]);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
}
