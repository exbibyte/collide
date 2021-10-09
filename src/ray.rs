use crate::bound::Bound;
use crate::bound_aabb::AxisAlignedBBox;
use crate::plane::*;
use crate::point::*;
use crate::ray_point_intersect;
use crate::ray_ray_intersect;
use crate::shape::*;
use crate::sphere::*;
use crate::vicinity::Vicinity;
use core::any::Any;
use lightmatrix::matrix::*;
use num_traits::{Float, NumAssign};

#[derive(Debug, Clone)]
pub struct Ray<T: NumAssign + Copy + Default + Float> {
    pub _ori: Matrix<T, 4, 1>,
    pub _dir: Matrix<T, 4, 1>,
    pub _bound: AxisAlignedBBox<T>,
    pub _vicinity: T,
}

impl<T> Ray<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn init(origin: &[T], dir: &[T]) -> Ray<T> {
        assert_eq!(origin.len(), 3);
        assert_eq!(dir.len(), 3);
        Ray {
            _ori: Matrix::from([[origin[0], origin[1], origin[2], T::one()]]).t(),
            _dir: Matrix::from([[dir[0], dir[1], dir[2], T::zero()]])
                .t()
                .normalize_l2(),
            _bound: AxisAlignedBBox::new(ShapeType::Ray, &[&origin[0..3], &dir[0..3]].concat()),
            _vicinity: T::from(1e-7).unwrap(),
        }
    }
}

impl<T> AnyBase for Ray<T>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<T> Shape<T> for Ray<T>
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
                let other_ray: &Ray<T> = match other.as_any().downcast_ref::<Self>() {
                    Some(b) => b,
                    None => {
                        panic!("cast to Ray failed");
                    }
                };
                ray_ray_intersect::intersect(self, other_ray)
            }
            ShapeType::Point => {
                let other_point: &Point<T> = match other.as_any().downcast_ref::<Point<T>>() {
                    Some(b) => b,
                    None => {
                        panic!("cast to Point failed");
                    }
                };
                ray_point_intersect::intersect(self, other_point)
            }
            ShapeType::Sphere => {
                let other_sphere: &Sphere<T> = match other.as_any().downcast_ref::<Sphere<T>>() {
                    Some(b) => b,
                    None => {
                        panic!("cast to Sphere failed");
                    }
                };

                let b_off = other_sphere._ori;
                let b_r = other_sphere._radius;

                let a_dir = self._dir;
                let a_off = self._ori;

                //sub in the ray equation into sphere equation
                // b := projection of relative offset onto ray direction
                // c := (minimal possible distance between sphere and ray origin )^2
                let relative_offset = a_off - b_off;
                let b = relative_offset.inner(&a_dir);
                let c = relative_offset.inner(&relative_offset) - b_r * b_r;

                if b > T::zero() && c > T::zero() {
                    //ray is outside of the sphere and points away from sphere
                    //thus no intersection occurs
                    return (false, None);
                }

                let d = b * b - c;
                if d < T::zero() {
                    //ray misses sphere
                    return (false, None);
                }

                let t1 = -b - d.sqrt();
                let t2 = -b + d.sqrt();

                let t = if t1 < T::zero() {
                    t2
                } else if t2 < T::zero() {
                    t1
                } else if t1 < t2 {
                    t1
                } else {
                    t2
                };

                return (true, Some((a_dir * t) + a_off));
            }
            ShapeType::Plane => {
                let other_plane: &Plane<T> = match other.as_any().downcast_ref::<Plane<T>>() {
                    Some(b) => b,
                    None => {
                        panic!("cast to Plane failed");
                    }
                };

                let b_off = other_plane._offset;
                let b_nor = other_plane._normal;

                //ray equation: r(t) = r.offset + r.dir * t
                //plane: p(x) = dot(normal, x-p.offset) = 0
                //p(x) = -dot(p.normal, p.offset) + dot(p.normal, x) = 0
                //substitution:
                // p(t) = -dot(p.fofset,p.normal) + dot(p.normal, r.offset + r.dir*t) = 0
                //      = -dot(p.fofset,p.normal) + dot(p.normal, r.offset) + t*dot(p.normal, r.dir) = 0
                //t = ( dot(p.offset, p.normal) - dot(p.normal, r.offset) )/ dot(p.normal, r.dir )
                let constant = b_off.inner(&b_nor);
                let numerator = constant - b_nor.inner(&self._ori);
                let denominator = b_nor.inner(&self._dir);
                if denominator == T::zero() {
                    //ray direction is colplaner to the plane
                    if constant == self._ori.inner(&b_nor) {
                        return (true, Some(self._ori.clone()));
                    } else {
                        return (false, None);
                    }
                } else if denominator > T::zero() {
                    //ray direction is not facing plane normal
                    return (false, None);
                }
                let t = numerator / denominator;
                if t < T::zero() {
                    return (false, None);
                }
                return (true, Some((self._dir * t) + self._ori));
            }
            _ => {
                unimplemented!();
            }
        }
    }
    fn get_support(&self, _v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>> {
        None
    }
}

impl<T> Vicinity<T> for Ray<T>
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
fn test_intersect_ray_ray() {
    //parallel rays, no intersection
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Ray::init(&[25f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);

        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for parallel rays, no intersection"),
        }
    }

    //colinear rays, intersection
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Ray::init(&[22f64, 2f64, 2f64], &[1f64, 1f64, 1f64]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == b._ori);
            }
            _ => panic!("unexpected result for parallel rays, no intersection"),
        }
    }

    //colinear rays, intersection
    {
        let a = Ray::init(&[25f64, 5f64, 5f64], &[1f64, 1f64, 1f64]);
        let b = Ray::init(&[22f64, 2f64, 2f64], &[1f64, 1f64, 1f64]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == a._ori);
            }
            _ => panic!("unexpected result for parallel rays, no intersection"),
        }
    }

    //rays, intersection
    {
        let a = Ray::init(&[5f64, 5f64, 0f64], &[-1f64, 0f64, 0f64]);
        let b = Ray::init(&[0f64, 0f64, 0f64], &[0f64, 1f64, 0f64]);

        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == Matrix::from([[0f64, 5f64, 0f64, 1f64]]).t());
            }
            _ => panic!("unexpected result for ray intersection"),
        }
    }

    //non-coplaner rays, no intersection
    {
        let a = Ray::init(&[5f64, 5f64, 2f64], &[-1f64, -1f64, 0f64]);
        let b = Ray::init(&[5f64, 5f64, 0f64], &[1f64, 1f64, 0f64]);

        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray intersection"),
        }
    }
}

#[test]
fn test_intersect_ray_point() {
    //ray point intersection
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Point::init(&[25f64, 5f64, 5f64]);
        match a.get_intersect(&b) {
            (true, Some(loc)) => {
                assert!(loc == b._ori);
            }
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
    //ray point no intersection, point behind ray origin and direction
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Point::init(&[15f64, -5f64, -5f64]);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray point no intersection, point behind ray"),
        }
    }
    //ray point no intersection
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Point::init(&[25f64, 5f64, 5.1f64]);
        match a.get_intersect(&b) {
            (false, None) => (),
            _ => panic!("unexpected result for ray point no intersection"),
        }
    }

    //repeat the above tests but invoking the method on point
    //ray point intersection
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Point::init(&[25f64, 5f64, 5f64]);
        match b.get_intersect(&a) {
            (true, Some(loc)) => {
                assert!(loc == b._ori);
            }
            _ => panic!("unexpected result for ray point intersection"),
        }
    }
    //ray point no intersection, point behind ray origin and direction
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Point::init(&[15f64, -5f64, -5f64]);
        match b.get_intersect(&a) {
            (false, None) => (),
            _ => panic!("unexpected result for ray point no intersection, point behind ray"),
        }
    }
    //ray point no intersection
    {
        let a = Ray::init(&[20f64, 0f64, 0f64], &[1f64, 1f64, 1f64]);
        let b = Point::init(&[25f64, 5f64, 5.1f64]);
        match b.get_intersect(&a) {
            (false, None) => (),
            _ => panic!("unexpected result for ray point no intersection"),
        }
    }
}
