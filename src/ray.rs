use crate::bound::Bound;
use crate::shape::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use num_traits::{Float, NumAssign};

use crate::bound_aabb::AxisAlignedBBox;
use lightmatrix::matrix::*;

use crate::ray_ray_intersect;

#[derive(Debug, Clone)]
pub struct Ray3<T: NumAssign + Copy + Default + Float> {
    pub _ori: Matrix<T, 4, 1>,
    pub _dir: Matrix<T, 4, 1>,
    pub _bound: AxisAlignedBBox<T>,
    pub _vicinity: T,
}

impl<T> Ray3<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn init(origin: &[T], dir: &[T]) -> Ray3<T> {
        assert_eq!(origin.len(), 3);
        assert_eq!(dir.len(), 3);
        Ray3 {
            _ori: Matrix::from([[origin[0], origin[1], origin[2], T::one()]]).t(),
            _dir: Matrix::from([[dir[0], dir[1], dir[2], T::zero()]])
                .t()
                .normalize_l2(),
            _bound: AxisAlignedBBox::new(ShapeType::Ray, &[&origin[0..3], &dir[0..3]].concat()),
            _vicinity: T::epsilon(),
        }
    }
}

impl<T> AnyBase for Ray3<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for Ray3<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn get_type(&self) -> ShapeType {
        ShapeType::Ray
    }
    fn get_bound(&self) -> &dyn Bound<T> {
        &self._bound
    }
    // this shall test for intersection of bounding shapes first before procedding to test intersection using algorithms of higher complexity
    fn get_intersect(&self, other: &dyn Shape<T>) -> (bool, Option<Matrix<T, 4, 1>>) {
        if !self.get_bound().intersect(other.get_bound()) {
            return (false, None);
        }
        match other.get_type() {
            ShapeType::Ray => {
                let other_ray: &Ray3<T> = match other.as_any().downcast_ref::<Self>() {
                    Some(b) => b,
                    None => {
                        panic!("cast to Ray3 failed");
                    }
                };
                ray_ray_intersect::intersect(self, other_ray)
            }
            // ShapeType::Point => {
            //     let other_shape_data = other.get_shape_data();
            //     let b_off = Matrix::from([[
            //         other_shape_data[0],
            //         other_shape_data[1],
            //         other_shape_data[2]]
            //     ]));
            //     let a_dir = &self._dir;
            //     let a_off = &self._ori;
            //     //a_dir * t + a_off = b_off
            //     //t = (b_off - a_off) / a_dir
            //     let t = &(&b_off - a_off) / a_dir;
            //     if !self.within_vicinity(t[0], t[1]) || !self.within_vicinity(t[1], t[2]) {
            //         return (false, None);
            //     } else {
            //         if t[0] >= 0f64 {
            //             return (true, Some(&(a_dir * t[0]) + a_off));
            //         } else {
            //             //the point is behind the ray origin and direction
            //             return (false, None);
            //         }
            //     }
            // }
            // ShapeType::Sphere => {
            //     let other_shape_data = other.get_shape_data();
            //     let ref b_off = Matrix::from([[
            //         other_shape_data[0],
            //         other_shape_data[1],
            //         other_shape_data[2]]
            //     ]));
            //     let b_r = other_shape_data[3];

            //     let ref a_dir = self._dir;
            //     let ref a_off = self._ori;

            //     //sub in the ray equation into sphere equation
            //     // b := projection of relative offset onto ray direction
            //     // c := (minimal possible distance between sphere and ray origin )^2
            //     let relative_offset = a_off - b_off;
            //     let b = relative_offset.inner(a_dir);
            //     let c = relative_offset.inner(&relative_offset) - b_r * b_r;

            //     if b > 0f64 && c > 0f64 {
            //         //ray is outside of the sphere and points away from sphere
            //         //thus no intersection occurs
            //         return (false, None);
            //     }

            //     let d = b * b - c;
            //     if d < 0f64 {
            //         //ray misses sphere
            //         return (false, None);
            //     }

            //     let t1 = -b - d.sqrt();
            //     let t2 = -b + d.sqrt();

            //     let t = if t1 < 0f64 {
            //         t2
            //     } else if t2 < 0f64 {
            //         t1
            //     } else if t1 < t2 {
            //         t1
            //     } else {
            //         t2
            //     };

            //     return (true, Some(&(a_dir * t) + a_off));
            // }
            // ShapeType::Plane => {
            //     let other_shape_data = other.get_shape_data();
            //     let b_off = Matrix::from([[
            //         other_shape_data[0],
            //         other_shape_data[1],
            //         other_shape_data[2]]
            //     ]));
            //     let b_nor = Matrix::from([[
            //         other_shape_data[3],
            //         other_shape_data[4],
            //         other_shape_data[5]]
            //     ]));
            //     //ray equation: r(t) = r.offset + r.dir * t
            //     //plane: p(x) = dot(normal, x-p.offset) = 0
            //     //p(x) = -dot(p.normal, p.offset) + dot(p.normal, x) = 0
            //     //substitution:
            //     // p(t) = -dot(p.fofset,p.normal) + dot(p.normal, r.offset + r.dir*t) = 0
            //     //      = -dot(p.fofset,p.normal) + dot(p.normal, r.offset) + t*dot(p.normal, r.dir) = 0
            //     //t = ( dot(p.offset, p.normal) - dot(p.normal, r.offset) )/ dot(p.normal, r.dir )
            //     let constant = b_off.inner(&b_nor);
            //     let numerator = constant - b_nor.inner(&self._ori);
            //     let denominator = b_nor.inner(&self._dir);
            //     if denominator == 0f64 {
            //         //ray direction is colplaner to the plane
            //         if constant == self._ori.inner(&b_nor) {
            //             return (true, Some(self._ori.clone()));
            //         } else {
            //             return (false, None);
            //         }
            //     } else if denominator > 0f64 {
            //         //ray direction is not facing plane normal
            //         return (false, None);
            //     }
            //     let t = numerator / denominator;
            //     if t < 0f64 {
            //         return (false, None);
            //     }
            //     return (true, Some(&(&self._dir * t) + &self._ori));
            // }
            _ => {
                unimplemented!();
            }
        }
    }
    fn get_support(&self, _v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>> {
        None
    }
}

impl<T> Vicinity<T> for Ray3<T>
where
    T: NumAssign + Copy + Default + Float,
{
    fn set_vicinity(&mut self, epsilon: T) {
        self._vicinity = epsilon.abs();
    }
    fn within_vicinity(&self, a: T, b: T) -> bool {
        if a + self._vicinity > b && a - self._vicinity < b {
            true
        } else {
            false
        }
    }
}

#[test]
fn test_intersect_ray_ray() {
    //parallel rays, no intersection
    {
        let a = Ray3::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Ray3::init(&[25f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);

        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for parallel rays, no intersection"),
        }
    }

    //colinear rays, intersection
    {
        let a = Ray3::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Ray3::init(&[22f64, 2f64, 2f64], &[1f64, 1f64, 1f64]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == b._ori);
            }
            _ => panic!("unexpected result for parallel rays, no intersection"),
        }
    }

    //colinear rays, intersection
    {
        let a = Ray3::init(&[25f64, 5f64, 5f64], &[1f64, 1f64, 1f64]);
        let b = Ray3::init(&[22f64, 2f64, 2f64], &[1f64, 1f64, 1f64]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == a._ori);
            }
            _ => panic!("unexpected result for parallel rays, no intersection"),
        }
    }

    //rays, intersection
    {
        let a = Ray3::init(&[5f64, 5f64, 0f64], &[-1f64, 0f64, 0f64]);
        let b = Ray3::init(&[0f64, 0f64, 0f64], &[0f64, 1f64, 0f64]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == Matrix::from([[0f64, 5f64, 0f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray intersection"),
        }
    }

    //non-coplaner rays, no intersection
    {
        let a = Ray3::init(&[5f64, 5f64, 2f64], &[-1f64, -1f64, 0f64]);
        let b = Ray3::init(&[5f64, 5f64, 0f64], &[1f64, 1f64, 0f64]);

        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray intersection"),
        }
    }
}
