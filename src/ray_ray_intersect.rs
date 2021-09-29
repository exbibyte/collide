use crate::ray::*;
use crate::vicinity::Vicinity;
use lightmatrix::matrix::*;
use num_traits::{Float, NumAssign};

pub(crate) fn intersect<T>(ray_a: &Ray<T>, ray_b: &Ray<T>) -> (bool, Option<Matrix<T, 4, 1>>)
where
    T: NumAssign + Copy + Default + Float,
{
    let a_dir = ray_a._dir;
    let b_dir = ray_b._dir;

    let a_off = ray_a._ori;
    let b_off = ray_b._ori;

    let c = b_dir - a_dir;
    let v = a_dir.cross(&b_dir);

    let dot_v_c = v.inner(&c);
    if !ray_a.within_vicinity(dot_v_c, T::zero()) {
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
        if !ray_a.within_vicinity(triangle_area, T::zero()) {
            //no overlap
            // println!( "parallel but non-overlapping lines" );
            (false, None)
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
                (true, Some(a_off))
            } else {
                //intersection at offset of ray b
                (true, Some(ray_a._dir * distance + ray_a._ori))
            }
        }
    } else {
        //solvable intersection exists
        let numerator = d.cross(&b_dir);
        let t = numerator.norm_l2() / v.norm_l2();
        if t < T::zero() {
            (false, None)
        } else {
            (true, Some(ray_a._dir * t + ray_a._ori))
        }
    }
}
