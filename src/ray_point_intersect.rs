use crate::point::*;
use crate::ray::*;
use crate::vicinity::*;
use lightmatrix::matrix::*;
use num_traits::{Float, NumAssign};

pub(crate) fn intersect<T>(a: &Ray<T>, b: &Point<T>) -> (bool, Option<Matrix<T, 4, 1>>)
where
    T: NumAssign + Copy + Default + Float,
{
    let b_off = b._ori;
    let a_dir = a._dir;
    let a_off = a._ori;

    //a_dir * t + a_off = b_off
    //t = (b_off - a_off) / a_dir
    let t = (b_off - a_off) / a_dir;
    if !a.within_vicinity(t[[0, 0]], t[[1, 0]]) || !a.within_vicinity(t[[1, 0]], t[[2, 0]]) {
        (false, None)
    } else if t[[0, 0]] >= T::zero() {
        (true, Some((a_dir * t[[0, 0]]) + a_off))
    } else {
        //the point is behind the ray origin and direction
        (false, None)
    }
}
