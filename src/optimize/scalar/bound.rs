use num_traits::{Float, FromPrimitive, Pow};
use std::fmt::Debug;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum BoundOptimizerError {
    #[error("Invalid ratio, the ratio must be greater than 1 and smaller than 2")]
    InvalidRatio,
    #[error("Invalid tolerance must be greater than 1")]
    InvalidTolerance,
}

#[allow(dead_code)]
/// An iterative bound based approach for minimization of a function.
///
/// # Arguments
/// * func: Function to minimize
/// * xi: The two values characterizing the bound
/// * ratio: The ratio around a bound to split around
/// * tolerance: The tolerance requirement to determine convergence
/// * max_iter: The maximum number of iterations to loop over.
///
/// # Example
/// ```
/// use mathslib::optimize::scalar::bound_optimizers::bound_minimize;
/// use mathslib::generals::Decimal;
///
/// fn case_1(x: f64) -> f64{x*x + 6.0*x + 3.0	}
///
/// fn main() {
/// 	assert_eq!(bound_minimize::<f64>(case_1, -9.0, 1.0, 1.5, 1e-4, 2000).unwrap().round_dp(4), -3.0000)
///}
/// ```
pub fn bound_minimize<T: Float + PartialOrd + Debug + FromPrimitive>(
    func: fn(T) -> T,
    mut x1: T,
    mut x2: T,
    ratio: T,
    tolerance: T,
    max_iter: u32,
) -> Result<T, BoundOptimizerError> {
    if tolerance < T::from_f64(0.0).unwrap() {
        return Err(BoundOptimizerError::InvalidTolerance);
    }

    if ratio < T::from_f64(1.0).unwrap() {
        return Err(BoundOptimizerError::InvalidRatio);
    }

    let mut c = x2 - (x2 - x1) / ratio;
    let mut d = x1 + (x2 - x1) / ratio;

    for _ in 0..max_iter {
        if (x2 - x1).abs() < tolerance {
            break;
        }
        if func(c) < func(d) {
            x2 = d;
        } else {
            x1 = c;
        }
        c = x2 - (x2 - x1) / ratio;
        d = x1 + (x2 - x1) / ratio;
    }

    Ok((x2 + x1) / T::from(2.0).unwrap())
}

#[allow(dead_code)]
/// An iterative bound based approach for minimization of a function using the golden ratio.
///
/// # Arguments
/// * func: Function to minimize
/// * xi: The two values characterizing the bound
/// * tolerance: The tolerance requirement to determine convergence
/// * max_iter: The maximum number of iterations to loop over.
///
/// # Example
/// ```
/// use mathslib::optimize::scalar::bound_optimizers::bound_gr_minimize;
/// use mathslib::generals::Decimal;
///
/// fn case_1(x: f64) -> f64{x*x + 6.0*x + 3.0	}
///
/// fn main() {
/// 	assert_eq!(bound_gr_minimize::<f64>(case_1, -9.0, 1.0, 1e-4, 2000).unwrap().round_dp(4), -3.0000)
///}
/// ```
pub fn bound_gr_minimize<T: Float + PartialOrd + Debug + FromPrimitive>(
    func: fn(T) -> T,
    x1: T,
    x2: T,
    tolerance: T,
    max_iter: u32,
) -> Result<T, BoundOptimizerError> {
    bound_minimize(
        func,
        x1,
        x2,
        T::from_f64((5.0.pow(0.5) + 1.0) / 2.0).unwrap(),
        tolerance,
        max_iter,
    )
}
