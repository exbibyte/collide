//! based on reference tutorial from http://www.dyn4j.org/2010/04/gjk-gilbert-johnson-keerthi/

use crate::shape::*;
use lightmatrix::matrix::*;
use num_traits::{Float, NumAssign};

fn support<T>(a: &dyn Shape<T>, b: &dyn Shape<T>, v: &Matrix<T, 4, 1>) -> Option<Matrix<T, 4, 1>>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    let p0 = match a.get_support(&v) {
        Some(o) => o,
        _ => return None,
    };
    let v_oppose = *v * -T::one();
    let p1 = match b.get_support(&v_oppose) {
        Some(o) => o,
        _ => return None,
    };
    let p10 = p0 - p1;
    Some(p10)
}

fn pass_minkowski_origin<T>(last_vert: &Matrix<T, 4, 1>, support: &Matrix<T, 4, 1>) -> bool
where
    T: NumAssign + Copy + Default + Float,
{
    // println!( "last vert dot product: {}", last_vert.dot( &support ).unwrap() );
    last_vert.inner(&support) > T::zero()
}

fn contains_minkowski_origin<T>(
    simplex: &mut Vec<Matrix<T, 4, 1>>,
    support: &mut Matrix<T, 4, 1>,
) -> bool
where
    T: NumAssign + Copy + Default + Float,
{
    let a = simplex.last().unwrap().clone();
    let ao = a * -T::one();
    if simplex.len() == 3 {
        //triangle case
        let b = simplex[1];
        let c = simplex[0];
        let ab = b - a;
        let ac = c - a;
        let ab_normal = ac.cross(&ab).cross(&ab);
        let ac_normal = ab.cross(&ac).cross(&ac);
        if ab_normal.inner(&ao) >= T::zero() {
            //remove c and set new direction to ab_normal
            let simplex_new = vec![simplex[1], simplex[2]];
            *simplex = simplex_new;
            *support = ab_normal.clone();
        } else if ac_normal.inner(&ao) >= T::zero() {
            //remove b and set new direction to ac_normal
            let simplex_new = vec![simplex[0], simplex[2]];
            *simplex = simplex_new.clone();
            *support = ac_normal.clone();
        } else {
            //minkowski origin is enclosed by the triangle
            return true;
        }
    } else {
        //line segment case
        //set direction towards minkowski origin
        let b = simplex[0];
        let ab = b - a;
        let ab_normal = ab.cross(&ao).cross(&ab);
        if ab_normal.norm_l2() == T::zero() {
            return true;
        } else {
            *support = ab_normal.clone();
        }
    }
    return false;
}

pub fn query_intersect<T>(a: &dyn Shape<T>, b: &dyn Shape<T>) -> Option<bool>
where
    T: NumAssign + Copy + Default + Float + 'static,
{
    match (a.get_type(), b.get_type()) {
        (ShapeType::Sphere, ShapeType::Sphere) => {}
        //todo
        // (ShapeType::POINT,ShapeType::BOX) => {},
        // (ShapeType::BOX,ShapeType::POINT) => {},
        _ => {
            panic!("unsupported shape type");
        }
    }
    //set initial minkowski vertex from an arbitrary support vector
    let mut d = Matrix::from([[-T::one(), T::zero(), T::zero(), T::zero()]]).t();

    let mut simplex = vec![];
    {
        let sup = support(a, b, &d).unwrap();
        simplex.push(sup);
    }

    d = d * -T::one();
    loop {
        // println!( "support vector: {:?}", d );
        {
            let sup = support(a, b, &d).unwrap();
            simplex.push(sup);
        }
        assert!(simplex.len() <= 3, "simplex vertices count unexpected");
        // println!( "simplex len: {}", simplex.len() );
        if !pass_minkowski_origin(simplex.last().unwrap(), &d) {
            // println!( "new vert not pass origin" );
            return Some(false);
        } else {
            if contains_minkowski_origin(&mut simplex, &mut d) {
                return Some(true);
            }
        }
    }
}

#[test]
fn test_intersect_gjk_shape_support() {
    use crate::sphere::*;
    for i in 0..10 {
        for j in 0..10 {
            for k in 0..10 {
                let a = Sphere::init(&[-5f64, 2.5f64, 15f64], 5.5f64);
                let v_x = 0.2f64 * (i as f64);
                let v_y = 0.2f64 * (j as f64);
                let v_z = 0.2f64 * (k as f64);

                let v = Matrix::from([[v_x, v_y, v_z, 0.]]).t();

                match a.get_support(&v) {
                    Some(o) => {
                        let l = (v_x * v_x + v_y * v_y + v_z * v_z).sqrt();
                        assert_eq!(
                            o,
                            Matrix::from([[
                                -5f64 + v_x / l * 5.5f64,
                                2.5f64 + v_y / l * 5.5f64,
                                15f64 + v_z / l * 5.5f64,
                                1.
                            ]])
                            .t()
                        );
                    }
                    _ => {
                        if i != 0 || j != 0 || k != 0 {
                            panic!("unexpected result");
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn test_intersect_gjk_query_intersect_positive_0() {
    use crate::sphere::*;
    let a = Sphere::init(&[0f64, 0f64, 0f64], 5f64);
    let b = Sphere::init(&[7f64, 0f64, 0f64], 2.1f64);
    let ret = query_intersect(&a, &b);
    assert!(ret.expect("gjk return unexpected"));
}
#[test]
fn test_intersect_gjk_query_intersect_positive_1() {
    use crate::sphere::*;
    let a = Sphere::init(&[0f64, 5f64, 0f64], 5f64);
    let b = Sphere::init(&[0f64, 0f64, 0f64], 2f64);
    let ret = query_intersect(&a, &b);
    assert!(ret.expect("gjk return unexpected"));
}
#[test]
fn test_intersect_gjk_query_intersect_positive_2() {
    use crate::sphere::*;
    let a = Sphere::init(&[0f64, 5f64, 0f64], 10f64);
    let b = Sphere::init(&[1f64, 1f64, 0f64], 2f64);
    let ret = query_intersect(&a, &b);
    assert!(ret.expect("gjk return unexpected"));
}
#[test]
fn test_intersect_gjk_query_intersect_positive_3() {
    use crate::sphere::*;
    let a = Sphere::init(&[0f64, -4.999f64, 0f64], 5f64);
    let b = Sphere::init(&[0f64, 5f64, 0f64], 5f64);
    let ret = query_intersect(&a, &b);
    assert!(ret.expect("gjk return unexpected"));
}
//todo
// {
//     let a = Point3::init( &[ -9.9, 9.9, 9.9 ] );
//     let b = RecBox::init( &[ 0.,0.,0. ], 10. );
//     let ret = intersect_gjk::query_intersect( &a, &b );
//     assert!( ret.expect("gjk return unexpected") );
// }

#[test]
fn test_intersect_gjk_query_intersect_negative_0() {
    use crate::sphere::*;
    let a = Sphere::init(&[0f64, 0f64, 0f64], 5f64);
    let b = Sphere::init(&[7f64, 0f64, 0f64], 1.99f64);
    let ret = query_intersect(&a, &b);
    assert!(!ret.expect("gjk return unexpected"));
}
#[test]
fn test_intersect_gjk_query_intersect_negative_1() {
    use crate::sphere::*;
    let a = Sphere::init(&[0f64, 5f64, 0f64], 5f64);
    let b = Sphere::init(&[0f64, 0f64, 10f64], 2f64);
    let ret = query_intersect(&a, &b);
    assert!(!ret.expect("gjk return unexpected"));
}
#[test]
fn test_intersect_gjk_query_intersect_negative_2() {
    use crate::sphere::*;
    let a = Sphere::init(&[0f64, -5f64, 0f64], 5f64);
    let b = Sphere::init(&[0f64, 5f64, 0f64], 5f64);
    let ret = query_intersect(&a, &b);
    assert!(!ret.expect("gjk return unexpected"));
}
