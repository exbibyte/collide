use crate::bound::Bound;
use crate::bound_aabb::AxisAlignedBBox;
use crate::ray::*;
use crate::shape::*;
use crate::vicinity::*;
use core::any::Any;
use lightmatrix::matrix::*;
use num_traits::{Float, NumAssign};

pub(crate) fn intersect<T>(a: &Ray3<T>, b: &Ray3<T>) -> (bool, Option<Matrix<T, 4, 1>>)
where
    T: NumAssign + Copy + Default + Float,
{
    let a_dir = a._dir;
    let b_dir = b._dir;

    let a_off = a._ori;
    let b_off = b._ori;

    let c = b_dir - a_dir;
    let v = a_dir.cross(&b_dir);

    let dot_v_c = v.inner(&c);
    if !a.within_vicinity(dot_v_c, T::zero()) {
        //they are not in the same place, so no intersection occurs
        return (false, None);
    }
    //test for colinearity
    let d = b_off - a_off;
    if v.inner(&v) < T::epsilon() {
        //lines are parallel
        //check triangle area formed by points on ray a and b
        let point1 = a_dir;
        let point2 = b_off - a_off;
        let triangle_area = point1.cross(&point2).norm_l2();
        // println!( "triangle area: {}", triangle_area );
        if !a.within_vicinity(triangle_area, T::zero()) {
            //no overlap
            // println!( "parallel but non-overlapping lines" );
            return (false, None);
        } else {
            //lines are colinear
            let direction = if d.inner(&a_dir) < T::zero() {
                -T::one()
            } else {
                T::one()
            };
            let distance = direction * d.norm_l2() / a_dir.norm_l2();
            // println!( "colinear lines, distance: {}", distance );
            if distance < T::zero() {
                //intersection at offset of ray a, so clamp t to 0
                return (true, Some(a_off.clone()));
            } else {
                //intersection at offset of ray b
                return (true, Some(a._dir * distance + a._ori));
            }
        }
    } else {
        //solvable intersection exists
        let numerator = d.cross(&b_dir);
        let t = numerator.norm_l2() / v.norm_l2();
        if t < T::zero() {
            return (false, None);
        } else {
            return (true, Some(a._dir * t + a._ori));
        }
    }
}
