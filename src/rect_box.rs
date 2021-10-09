use crate::bound::Bound;
use crate::point::*;
use crate::shape::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use num_traits::{Float, NumAssign};

use crate::bound_aabb::AxisAlignedBBox;
use lightmatrix::matrix::*;

#[derive(Debug, Clone)]
pub struct RectBox<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub _ori: Matrix<T, 4, 1>,
    pub _size: T,
    pub _bound: AxisAlignedBBox<T>,
    pub _vicinity: T,
}

impl<T> RectBox<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn init(origin: &[T], size: T) -> RectBox<T> {
        assert!(origin.len() == 3);
        RectBox {
            _ori: Matrix::from([[origin[0], origin[1], origin[2], T::one()]]).t(),
            _size: size, //half of the length of box edge
            _bound: AxisAlignedBBox::new(ShapeType::Box, &[&origin[0..3], &[size]].concat()),
            _vicinity: T::epsilon(),
        }
    }
}

impl<T> AnyBase for RectBox<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for RectBox<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType {
        ShapeType::Box
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
                ShapeType::Point => {
                    //covered by bbox test
                    let other_point: &Point<T> = match other.as_any().downcast_ref::<Point<T>>() {
                        Some(b) => b,
                        None => {
                            panic!("cast to Point failed");
                        }
                    };
                    let b_off = other_point._ori;
                    return (true, Some(b_off));
                }
                _ => {
                    unimplemented!();
                }
            }
        }
    }
    fn get_support(&self, v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>> {
        if v.norm_l2() != T::zero() {
            //get a furthest point in the given direction v
            let points = [
                Matrix::from([[self._size, self._size, self._size, T::one()]]).t(),
                Matrix::from([[-self._size, self._size, self._size, T::one()]]).t(),
                Matrix::from([[self._size, -self._size, self._size, T::one()]]).t(),
                Matrix::from([[-self._size, -self._size, self._size, T::one()]]).t(),
                Matrix::from([[self._size, self._size, -self._size, T::one()]]).t(),
                Matrix::from([[-self._size, self._size, -self._size, T::one()]]).t(),
                Matrix::from([[self._size, -self._size, -self._size, T::one()]]).t(),
                Matrix::from([[-self._size, -self._size, -self._size, T::one()]]).t(),
            ];

            let furthest = points
                .iter()
                .map(|x| x.inner(v))
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .unwrap();

            let o = self._ori + points[furthest.0];
            Some(o)
        } else {
            None
        }
    }
}

impl<T> Vicinity<T> for RectBox<T>
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
fn test_intersect_rectbox_point_0() {
    //intersection
    {
        let a = Point::init(&[-9.9, 9.9, 9.9]);
        let b = RectBox::init(&[0., 0., 0.], 10.);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                println!("{:?}", loc);
                println!("{:?}", a._ori);
                assert_eq!(loc, a._ori);
                //     assert!(loc.is_equal(&a._ori, 0.0001f64).unwrap());
            }
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
}
#[test]
fn test_intersect_rectbox_point_1() {
    //no intersection
    {
        let a = Point::init(&[-9.9, 9.9, -10.1]);
        let b = RectBox::init(&[0., 0., 0.], 10.);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
}
