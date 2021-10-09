use crate::bound::Bound;
use crate::bound_aabb::AxisAlignedBBox;
use crate::point::*;
#[cfg(test)]
use crate::ray::*;
use crate::shape::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use lightmatrix::matrix::*;
use num_traits::{Float, NumAssign};

#[derive(Debug, Clone)]
pub struct Sphere<T: NumAssign + Copy + Default + Float> {
    pub _ori: Matrix<T, 4, 1>,
    pub _radius: T,
    pub _bound: AxisAlignedBBox<T>,
    pub _vicinity: T,
}

impl<T> Sphere<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn init(origin: &[T], r: T) -> Sphere<T> {
        assert!(origin.len() == 3);
        Sphere {
            _ori: Matrix::from([[origin[0], origin[1], origin[2], T::one()]]).t(),
            _radius: r,
            _bound: AxisAlignedBBox::new(ShapeType::Sphere, &[&origin[0..3], &[r]].concat()),
            _vicinity: T::from(1e-7).unwrap(),
        }
    }
}

impl<T> AnyBase for Sphere<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for Sphere<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType {
        ShapeType::Sphere
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
                ShapeType::Sphere => {
                    let other_sphere: &Sphere<T> = match other.as_any().downcast_ref::<Self>() {
                        Some(b) => b,
                        None => {
                            panic!("cast to Sphere failed");
                        }
                    };

                    let b_off = other_sphere._ori;
                    let a_r = self._radius;
                    let b_r = other_sphere._radius;

                    let a_off = self._ori;
                    let c = b_off - a_off;
                    let d = c.norm_l2();
                    if d > b_r + a_r {
                        return (false, None);
                    } else {
                        //calculate a mid point average
                        let f = a_r / (a_r + b_r);
                        let g = c * f;
                        return (true, Some(a_off + g));
                    }
                }
                ShapeType::Ray => {
                    //see Ray for ray sphere intersection
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
                    let d = b_off - self._ori;
                    for i in 0..3 {
                        if d[[i, 0]] > self._radius {
                            return (false, None);
                        }
                    }
                    return (true, Some(b_off));
                }
                // ShapeType::Plane => {
                //     let other_shape_data = other.get_shape_data();
                //     let b_off = Mat3x1 {
                //         _val: [
                //             other_shape_data[0],
                //             other_shape_data[1],
                //             other_shape_data[2],
                //         ],
                //     };
                //     let b_nor = Mat3x1 {
                //         _val: [
                //             other_shape_data[3],
                //             other_shape_data[4],
                //             other_shape_data[5],
                //         ],
                //     };
                //     //x = -plane_normal * t + sphere_center
                //     //dot( plane_normal, x ) = dot( plane_normal, plane_offset ) = k
                //     //substitution:
                //     //dot( plane_normal, -plane_normal * t + sphere_center ) = k
                //     //-t + dot( plane_normal, sphere_center ) = k
                //     //t = dot( plane_normal, sphere_center ) - k

                //     let k = b_nor.dot(&b_off).unwrap();
                //     let t = b_nor.dot(&self._ori).unwrap() - k;
                //     if t > self._radius {
                //         return (false, None);
                //     } else {
                //         return (
                //             true,
                //             Some(b_nor.scale(-t).unwrap().plus(&self._ori).unwrap()),
                //         );
                //     }
                // }
                _ => {
                    unimplemented!();
                }
            }
        }
    }
    fn get_support(&self, v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>> {
        if v.norm_l2() != T::zero() {
            let v_adjusted = v.normalize_l2() * self._radius;
            let o = self._ori + v_adjusted;
            Some(o)
        } else {
            None
        }
    }
}

impl<T> Vicinity<T> for Sphere<T>
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
fn test_intersect_sphere_sphere_0() {
    //Sphere Sphere intersection
    {
        let a = Sphere::init(&[10f64, 0f64, 0f64], 5f64);
        let b = Sphere::init(&[20f64, 0f64, 0f64], 5f64);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, Matrix::from([[15f64, 0f64, 0f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray sphere intersection"),
        }
    }
}
#[test]
fn test_intersect_sphere_sphere_1() {
    //Sphere Sphere intersection
    {
        let a = Sphere::init(&[10f64, 0f64, 0f64], 5f64);
        let b = Sphere::init(&[13f64, 4f64, 0f64], 5f64);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                println!("loc: {:?}", loc);
                assert_eq!(loc, Matrix::from([[11.5f64, 2f64, 0f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray sphere intersection"),
        }
    }
}
#[test]
fn test_intersect_sphere_sphere_2() {
    //Sphere Sphere no intersection
    {
        let a = Sphere::init(&[10f64, 0f64, 0f64], 5f64);
        let b = Sphere::init(&[20f64, 0.1f64, 0f64], 5f64);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray sphere intersection"),
        }
    }
}

#[test]
fn test_intersect_sphere_point_0() {
    //sphere point intersection
    {
        let a = Sphere::init(&[10f64, 0f64, 0f64], 5f64);
        let b = Point::init(&[8f64, 2f64, 3f64]);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, Matrix::from([[8f64, 2f64, 3f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
}
#[test]
fn test_intersect_sphere_point_1() {
    //sphere point intersection
    {
        let a = Sphere::init(&[10f64, 0f64, 0f64], 5f64);
        let b = Point::init(&[10f64, 5f64, 0f64]);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, Matrix::from([[10f64, 5f64, 0f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
}
#[test]
fn test_intersect_sphere_point_2() {
    //sphere point no intersection
    {
        let a = Sphere::init(&[10f64, 0f64, 0f64], 5f64);
        let b = Point::init(&[0f64, 5.1f64, 0f64]);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray point no intersection"),
        }
    }
}

#[test]
fn test_intersect_ray_sphere_0() {
    //Ray Sphere intersection
    {
        let a = Ray::init(&[5f64, 0f64, 0f64], &[1f64, 0f64, 0f64]);
        let b = Sphere::init(&[20f64, 0f64, 0f64], 5f64);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, Matrix::from([[15f64, 0f64, 0f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray sphere intersection"),
        }
    }
}
#[test]
fn test_intersect_ray_sphere_1() {
    //Ray Sphere no intersection, opposing direction
    {
        let a = Ray::init(&[5f64, 0f64, 0f64], &[-1f64, 0f64, 0f64]);
        let b = Sphere::init(&[20f64, 0f64, 0f64], 5f64);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray sphere no intersection"),
        }
    }
}
#[test]
fn test_intersect_ray_sphere_2() {
    //Ray Sphere intersection, at edge
    {
        let a = Ray::init(&[30f64, 10f64, 10f64], &[-1f64, 0f64, 0f64]);
        let b = Sphere::init(&[20f64, 10f64, 10f64], 5f64);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert_eq!(loc, Matrix::from([[25f64, 10f64, 10f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray sphere intersection"),
        }
    }
}
#[test]
fn test_intersect_ray_sphere_3() {
    //Ray Sphere intersection, oblique angle
    {
        let a = Ray::init(&[30f64, 10f64, 10f64], &[-1f64, -1f64, -1f64]);
        let b = Sphere::init(&[20f64, 0f64, 0f64], 5f64);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                let c = (25f64 / 3f64).sqrt();
                println!("loc: {:?}", loc);

                assert_matrix_approx_eq_float(
                    &loc,
                    &Matrix::from([[20f64 + c, 0f64 + c, 0f64 + c, 1f64]]).t(),
                    1e-7,
                );
            }
            _ => panic!("unexpected result for ray sphere intersection"),
        }
    }
}
#[test]
fn test_intersect_ray_sphere_4() {
    //Ray Sphere no intersection
    {
        let a = Ray::init(&[30f64, 10f64, 10f64], &[-1f64, 0f64, -1f64]);
        let b = Sphere::init(&[20f64, 0f64, 0f64], 5f64);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray sphere no intersection"),
        }
    }
}
