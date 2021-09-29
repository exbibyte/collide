use num_traits::{Float, NumAssign};
use std::cmp;

use crate::bound::Bound;
use crate::bound::BoundType;
use crate::shape::ShapeType;

#[derive(Debug, Clone, Copy)]
pub struct AxisAlignedBBox<T: NumAssign + Copy + Default + Float> {
    pub bound_lower: [T; 3],
    pub bound_upper: [T; 3],
}

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl<T> AxisAlignedBBox<T>
where
    T: NumAssign + Copy + Default + Float,
{
    pub fn new(shape_type: ShapeType, vals: &[T]) -> AxisAlignedBBox<T> {
        match shape_type {
            ShapeType::Ray => {
                assert!(vals.len() == 6);
                let mut bounds = [(T::zero(), T::zero()); 3];
                for i in 0..3 {
                    let b = if vals[3 + i] > T::zero() {
                        (vals[i], T::infinity())
                    } else if vals[3 + i] < T::zero() {
                        (T::neg_infinity(), vals[i])
                    } else {
                        (vals[i], vals[i])
                    };
                    bounds[i] = b;
                }
                AxisAlignedBBox {
                    bound_lower: [bounds[0].0, bounds[1].0, bounds[2].0],
                    bound_upper: [bounds[0].1, bounds[1].1, bounds[2].1],
                }
            }
            ShapeType::Point => {
                assert!(vals.len() == 3);
                AxisAlignedBBox {
                    bound_lower: [vals[0], vals[1], vals[2]],
                    bound_upper: [vals[0], vals[1], vals[2]],
                }
            }
            ShapeType::Sphere => {
                assert!(vals.len() == 4);
                AxisAlignedBBox {
                    bound_lower: [vals[0] - vals[3], vals[1] - vals[3], vals[2] - vals[3]],
                    bound_upper: [vals[0] + vals[3], vals[1] + vals[3], vals[2] + vals[3]],
                }
            }
            ShapeType::Plane => {
                assert!(vals.len() == 6);
                AxisAlignedBBox {
                    bound_lower: [T::neg_infinity(); 3],
                    bound_upper: [T::infinity(); 3],
                }
            }
            ShapeType::Box => {
                assert!(vals.len() == 4);
                AxisAlignedBBox {
                    bound_lower: [vals[0] - vals[3], vals[1] - vals[3], vals[2] - vals[3]],
                    bound_upper: [vals[0] + vals[3], vals[1] + vals[3], vals[2] + vals[3]],
                }
            }
            ShapeType::Rect => {
                assert!(vals.len() == 6);
                AxisAlignedBBox {
                    bound_lower: [vals[0], vals[1], vals[2]],
                    bound_upper: [vals[3], vals[4], vals[5]],
                }
            }
            ShapeType::Frustum => {
                unimplemented!();
            }
            _ => {
                unimplemented!();
            }
        }
    }
    pub fn get_longest_axis(&self) -> (Axis, T) {
        let dx = (Axis::X, self.bound_upper[0] - self.bound_lower[0]);
        let dy = (Axis::Y, self.bound_upper[1] - self.bound_lower[1]);
        let dz = (Axis::Z, self.bound_upper[2] - self.bound_lower[2]);
        let longest = [dx, dy, dz]
            .iter()
            .cloned()
            .max_by(|x, y| {
                if x.1 < y.1 {
                    cmp::Ordering::Less
                } else if x.1 < y.1 {
                    cmp::Ordering::Greater
                } else {
                    cmp::Ordering::Equal
                }
            })
            .unwrap();
        longest
    }
}

impl<T> Bound<T> for AxisAlignedBBox<T>
where
    T: NumAssign + Copy + Default + Float,
{
    fn get_type(&self) -> BoundType {
        BoundType::AxisAlignBox
    }
    fn intersect(&self, other: &dyn Bound<T>) -> bool {
        match other.get_type() {
            BoundType::AxisAlignBox => {
                let a_bounds = self.get_bound_data();
                let b_bounds = other.get_bound_data();

                let a_lower = &a_bounds[0..3];
                let a_upper = &a_bounds[3..6];
                let b_lower = &b_bounds[0..3];
                let b_upper = &b_bounds[3..6];

                for i in 0..3 {
                    if a_lower[i] > b_upper[i] || a_upper[i] < b_lower[i] {
                        return false;
                    }
                }
                return true;
            }
            _ => {
                unimplemented!();
            }
        }
    }
    fn get_shortest_separation(&self, _other: &dyn Bound<T>) -> T {
        unimplemented!();
    }
    fn get_bound_data(&self) -> [T; 32] {
        let mut arr = [T::zero(); 32];
        for i in 0..3 {
            arr[i] = self.bound_lower[i];
        }
        for i in 0..3 {
            arr[i + 3] = self.bound_upper[i];
        }
        arr
    }
    fn get_union(&mut self, bounds: &[&dyn Bound<T>]) {
        self.bound_lower = [T::infinity(); 3];
        self.bound_upper = [T::neg_infinity(); 3];
        for i in bounds {
            match i.get_type() {
                BoundType::AxisAlignBox => (),
                _ => {
                    unimplemented!();
                }
            }
            let b = i.get_bound_data();
            let b_lower = &b[0..3];
            let b_upper = &b[3..6];
            for j in 0..3 {
                self.bound_lower[j] = self.bound_lower[j].min(b_lower[j]);
                self.bound_upper[j] = self.bound_upper[j].max(b_upper[j]);
            }
        }
    }
    fn get_centroid(&self) -> [T; 3] {
        match self.get_type() {
            BoundType::AxisAlignBox => {
                let b = self.get_bound_data();
                let b_lower = &b[0..3];
                let b_upper = &b[3..6];
                return [
                    (b_lower[0] + b_upper[0]) / T::from(2.).unwrap(),
                    (b_lower[1] + b_upper[1]) / T::from(2.).unwrap(),
                    (b_lower[2] + b_upper[2]) / T::from(2.).unwrap(),
                ];
            }
            _ => {
                unimplemented!();
            }
        }
    }
}

impl<T> Default for AxisAlignedBBox<T>
where
    T: NumAssign + Copy + Default + Float,
{
    fn default() -> AxisAlignedBBox<T> {
        AxisAlignedBBox {
            bound_lower: [T::neg_infinity(); 3],
            bound_upper: [T::infinity(); 3],
        }
    }
}
