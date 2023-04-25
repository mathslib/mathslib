use std::fmt::Debug;
use std::ops::AddAssign;
use num_traits::{abs, Float, FromPrimitive};
use crate::generals::differential_methods::{central_finite_difference, FiniteDifferenceError};
use thiserror::Error;


#[derive(Error, Debug)]
pub enum NewtonRaphsonError {
	#[error("Failed to calculate the derivative")]
	DerivativeError {
		#[from] source: FiniteDifferenceError
	},
	#[error("InvalidTolerance")]
	InvalidTolerance
}



fn newton_raphson<T: Float + Debug + AddAssign + FromPrimitive + num_traits::Signed>(func: fn(T) -> T, x0: T, tolerance: T, max_iter: u32, h: T) -> Result<T, NewtonRaphsonError> {
	// Validate the tolerance
	if tolerance < T::from_f64(0.0).unwrap() {
		return Err(NewtonRaphsonError::InvalidTolerance)
	}

	let mut x0= x0;
	let mut old_val: T = x0 + tolerance;

	for _ in 0..max_iter {

		// TODO: Create abstraction for higher order derivatives
		// FIXME: The following unwrap must be handled in the above abstraction
		x0 = x0 - central_finite_difference(func, x0, h, 1, 1)? / central_finite_difference(func, x0, h, 1, 2)?;

		if abs(old_val - x0) < tolerance {
			break
		}

		old_val = x0
	}

	Ok(x0)
}

#[cfg(test)]
mod test{
	use num_traits::Pow;
	use crate::optimize::scalar::newtonraphson::newton_raphson;

	fn case1(x: f64) -> f64 {
		x.pow(2) + 6.0*x + 3.0
	}

	#[test]
	fn test_newtonraphson() {
		assert_eq!(newton_raphson(case1, 0.0, 1e-5, 100, 1e5).unwrap(), -3.0)
	}
}

