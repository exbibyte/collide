use crate::bound::Bound;
use crate::line::*;
use crate::plane::*;
use crate::point::*;
use crate::ray::*;
use crate::shape::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use num_traits::{Float, NumAssign};

use crate::bound_aabb::AxisAlignedBBox;
use lightmatrix::matrix::*;

#[derive(Debug, Clone, Default)]
pub struct TriPrism<T>
where
    T: NumAssign + Copy + Default + Float,
{
    ///base
    pub _tri_base: [Matrix<T, 4, 1>; 3],

    ///base + height offset in normal direction of base
    pub _tri_base2: [Matrix<T, 4, 1>; 3],

    ///normal of the triangle base, scaled with height
    pub _normal_height: Matrix<T, 4, 1>,

    pub _bound: AxisAlignedBBox<T>,

    pub _vicinity: T,
}

impl<T> TriPrism<T>
where
    T: NumAssign + Copy + Default + Float,
{
    /// initialize with tribase: base vertices in ccw order
    pub fn init(tri_base: &[T], height: T) -> TriPrism<T> {
        assert!(tri_base.len() == 9);

        let v0 = Matrix::from([[tri_base[0], tri_base[1], tri_base[2], T::one()]]).t();

        let v1 = Matrix::from([[tri_base[3], tri_base[4], tri_base[5], T::one()]]).t();

        let v2 = Matrix::from([[tri_base[6], tri_base[7], tri_base[8], T::one()]]).t();

        let d1 = v1 - v0;
        let d2 = v2 - v0;
        let normal = d1.cross(&d2).normalize_l2();
        let h_offset = normal * height;

        let v00 = v0 + h_offset;
        let v11 = v1 + h_offset;
        let v22 = v2 + h_offset;

        let base = [v0, v1, v2];
        let base2 = [v00, v11, v22];

        use std::cmp::Ordering::*;

        let xs = [
            base[0][[0, 0]],
            base[1][[0, 0]],
            base[2][[0, 0]],
            base2[0][[0, 0]],
            base2[1][[0, 0]],
            base2[2][[0, 0]],
        ];

        let ys = [
            base[0][[1, 0]],
            base[1][[1, 0]],
            base[2][[1, 0]],
            base2[0][[1, 0]],
            base2[1][[1, 0]],
            base2[2][[1, 0]],
        ];

        let zs = [
            base[0][[2, 0]],
            base[1][[2, 0]],
            base[2][[2, 0]],
            base2[0][[2, 0]],
            base2[1][[2, 0]],
            base2[2][[2, 0]],
        ];

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

        TriPrism {
            _tri_base: base,
            _tri_base2: base2,
            _normal_height: h_offset,
            _bound: AxisAlignedBBox::new(
                ShapeType::Rect,
                &[x_min, y_min, z_min, x_max, y_max, z_max],
            ),
            _vicinity: T::epsilon(),
        }
    }
}

impl<T> AnyBase for TriPrism<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for TriPrism<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType {
        ShapeType::TriPrism
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
                    let other_point: &Point<T> = match other.as_any().downcast_ref::<Point<T>>() {
                        Some(b) => b,
                        None => {
                            panic!("cast to Point failed");
                        }
                    };

                    let other_point = other_point._ori;

                    //test point aginst 5 half spaces from facets of the tri_prism to determine if point is inside the tri_prism

                    let n = self._normal_height;

                    let tests = vec![
                        (self._tri_base[0], n * -T::one()),
                        (self._tri_base2[0], n),
                        (
                            self._tri_base[0],
                            (self._tri_base[1] - self._tri_base[0]).cross(&n),
                        ),
                        (
                            self._tri_base[1],
                            (self._tri_base[2] - self._tri_base[1]).cross(&n),
                        ),
                        (
                            self._tri_base[2],
                            (self._tri_base[0] - self._tri_base[2]).cross(&n),
                        ),
                    ];

                    let is_inside = tests
                        .iter()
                        .all(|(vert, normal)| !((other_point - *vert).inner(normal) > T::zero()));

                    if is_inside {
                        (true, Some(other_point))
                    } else {
                        (false, None)
                    }
                }
                ShapeType::Line => {
                    let other_point: &Line<T> = match other.as_any().downcast_ref::<Line<T>>() {
                        Some(b) => b,
                        None => {
                            panic!("cast to Point failed");
                        }
                    };

                    let a = other_point._a;
                    let b = other_point._b;

                    //test points aginst 5 half spaces from facets of the tri_prism to determine if point is inside the tri_prism

                    let n = self._normal_height;

                    let tests = vec![
                        (self._tri_base[0], n * -T::one()),
                        (self._tri_base2[0], n),
                        (
                            self._tri_base[0],
                            (self._tri_base[1] - self._tri_base[0]).cross(&n),
                        ),
                        (
                            self._tri_base[1],
                            (self._tri_base[2] - self._tri_base[1]).cross(&n),
                        ),
                        (
                            self._tri_base[2],
                            (self._tri_base[0] - self._tri_base[2]).cross(&n),
                        ),
                    ];

                    let a_is_inside = tests
                        .iter()
                        .all(|(vert, normal)| !((a - *vert).inner(normal) > T::zero()));
                    let b_is_inside = tests
                        .iter()
                        .all(|(vert, normal)| !((b - *vert).inner(normal) > T::zero()));

                    if a_is_inside {
                        return (true, Some(a));
                    } else if b_is_inside {
                        return (true, Some(b));
                    }

                    //continue test using ray plane intersection

                    let v = b - a;
                    let mag = v.norm_l2();

                    let r = Ray::init(
                        &[a[[0, 0]], a[[1, 0]], a[[2, 0]]],
                        &[v[[0, 0]], v[[1, 0]], v[[2, 0]]],
                    );

                    let n1 = self._normal_height;
                    let n0 = n1 * -T::one();
                    let n2 = (self._tri_base[1] - self._tri_base[0]).cross(&n1);
                    let n3 = (self._tri_base[2] - self._tri_base[1]).cross(&n1);
                    let n4 = (self._tri_base[0] - self._tri_base[2]).cross(&n1);

                    let facets = vec![
                        Plane::init(
                            &[
                                self._tri_base[0][[0, 0]] as T,
                                self._tri_base[0][[1, 0]] as T,
                                self._tri_base[0][[2, 0]] as T,
                            ],
                            &[n0[[0, 0]] as T, n0[[1, 0]] as T, n0[[2, 0]] as T],
                        ),
                        Plane::init(
                            &[
                                self._tri_base2[0][[0, 0]] as T,
                                self._tri_base2[0][[1, 0]] as T,
                                self._tri_base2[0][[2, 0]] as T,
                            ],
                            &[n1[[0, 0]] as T, n1[[1, 0]] as T, n1[[2, 0]] as T],
                        ),
                        Plane::init(
                            &[
                                self._tri_base[0][[0, 0]] as T,
                                self._tri_base[0][[1, 0]] as T,
                                self._tri_base[0][[2, 0]] as T,
                            ],
                            &[n2[[0, 0]] as T, n2[[1, 0]] as T, n2[[2, 0]] as T],
                        ),
                        Plane::init(
                            &[
                                self._tri_base[1][[0, 0]] as T,
                                self._tri_base[1][[1, 0]] as T,
                                self._tri_base[1][[2, 0]] as T,
                            ],
                            &[n3[[0, 0]] as T, n3[[1, 0]] as T, n3[[2, 0]] as T],
                        ),
                        Plane::init(
                            &[
                                self._tri_base[2][[0, 0]] as T,
                                self._tri_base[2][[1, 0]] as T,
                                self._tri_base[2][[2, 0]] as T,
                            ],
                            &[n4[[0, 0]] as T, n4[[1, 0]] as T, n4[[2, 0]] as T],
                        ),
                    ];

                    let mut intersect_point = None;
                    let mut is_inside = false;
                    for i in facets.iter() {
                        let res = r.get_intersect(i);
                        if res.0 {
                            let collide_point = res.1.unwrap();
                            let mag2 = (collide_point - a).norm_l2();

                            //one more check necesary for the candidate collision point
                            let is_point_inside = tests.iter().all(|(vert, normal)| {
                                !((collide_point - *vert).inner(normal) > T::zero())
                            });

                            if !is_point_inside || mag2 > mag {
                                continue;
                            } else {
                                is_inside = true;
                                intersect_point = res.1;
                                break;
                            }
                        }
                    }

                    if is_inside {
                        (true, intersect_point)
                    } else {
                        (false, None)
                    }
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
                self._tri_base[0],
                self._tri_base[1],
                self._tri_base[2],
                self._tri_base2[0],
                self._tri_base2[1],
                self._tri_base2[2],
            ];

            let furthest = points
                .iter()
                .map(|x| x.inner(v))
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .unwrap();

            let o = points[furthest.0].clone();

            Some(o)
        } else {
            None
        }
    }
}

impl<T> Vicinity<T> for TriPrism<T>
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
fn test_intersect_triprism_line_0() {
    //intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0.25, 0.25, 0.], &[1., 1., 0.]);

        match a.get_intersect(&b) {
            (true, Some(_loc)) => {}
            _ => panic!("unexpected result for triprism line intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_line_1() {
    //intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0.5, -50., 0.], &[0.5, 50., 0.]);

        match a.get_intersect(&b) {
            (true, Some(_loc)) => {}
            _ => panic!("unexpected result for triprism line intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_line_2() {
    //intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0.25, 0.25, 0.5], &[0.26, 0.26, 0.5]);

        match a.get_intersect(&b) {
            (true, Some(_loc)) => {}
            _ => panic!("unexpected result for triprism line intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_line_3() {
    //intersection, flipped
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0.25, 0.25, 0.], &[1., 1., 0.]);

        match b.get_intersect(&a) {
            (true, Some(_loc)) => {}
            _ => panic!("unexpected result for triprism line intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_line_4() {
    //intersection, flipped
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0.5, -50., 0.], &[0.5, 50., 0.]);

        match b.get_intersect(&a) {
            (true, Some(_loc)) => {}
            _ => panic!("unexpected result for triprism line intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_line_5() {
    //intersection, flipped
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0.25, 0.25, 0.5], &[0.26, 0.26, 0.5]);

        match b.get_intersect(&a) {
            (true, Some(_loc)) => {}
            _ => panic!("unexpected result for triprism line intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_line_6() {
    //no intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0.25, 0.25, 1.5], &[1., 1., 1.5]);

        match a.get_intersect(&b) {
            (true, Some(_loc)) => {
                panic!("unexpected result for triprism line intersection");
            }
            _ => {}
        }
    }
}
#[test]
fn test_intersect_triprism_line_7() {
    //no intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Line::init(&[0., -5., 0.5], &[50., 45., 0.5]);

        match a.get_intersect(&b) {
            (true, Some(_loc)) => {
                panic!("unexpected result for triprism line intersection");
            }
            _ => {}
        }
    }
}

#[test]
fn test_intersect_triprism_point_0() {
    //intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.25, 0.25, 0.]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, b._ori);
                // assert!(loc.is_equal(&b._ori, 0.0001f64).unwrap());
            }
            _ => panic!("unexpected result for triprism point intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_point_1() {
    //intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.25, 0.25, 0.5]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, b._ori);
                // assert!(loc.is_equal(&b._ori, 0.0001f64).unwrap());
            }
            _ => panic!("unexpected result for triprism point intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_point_2() {
    //intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.25, 0.25, 1.]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, b._ori);
                // assert!(loc.is_equal(&b._ori, 0.0001f64).unwrap());
            }
            _ => panic!("unexpected result for triprism point intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_point_3() {
    //intersection, flipped
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.25, 0.25, 0.]);

        match b.get_intersect(&a) {
            (true, Some(loc)) => {
                assert_eq!(loc, b._ori);
                // assert!(loc.is_equal(&b._ori, 0.0001f64).unwrap());
            }
            _ => panic!("unexpected result for triprism point intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_point_4() {
    //intersection, flipped
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.25, 0.25, 0.5]);

        match b.get_intersect(&a) {
            (true, Some(loc)) => {
                assert_eq!(loc, b._ori);
                // assert!(loc.is_equal(&b._ori, 0.0001f64).unwrap());
            }
            _ => panic!("unexpected result for triprism point intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_point_5() {
    //intersection, flipped
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.25, 0.25, 1.]);

        match b.get_intersect(&a) {
            (true, Some(loc)) => {
                assert_eq!(loc, b._ori);
                // assert!(loc.is_equal(&b._ori, 0.0001f64).unwrap());
            }
            _ => panic!("unexpected result for triprism point intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_point_6() {
    //no intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.25, 0.25, 1.001]);

        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for triprism point no intersection"),
        }
    }
}
#[test]
fn test_intersect_triprism_point_7() {
    //no intersection
    {
        let a = TriPrism::init(&[0., 0., 0., 1., 0., 0., 1., 1., 0.], 1.);

        let b = Point::init(&[0.5, 0.55, 0.5]);

        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for triprism point no intersection"),
        }
    }
}
