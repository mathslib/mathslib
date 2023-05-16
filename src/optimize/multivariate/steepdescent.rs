use std::fmt::Debug;
use std::ops::AddAssign;
use num_traits::{Float, FromPrimitive};
use crate::generals::differential_methods::{multivariate_central_finite_difference, FiniteDifferenceError};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SteepDescentError {
	#[error("Failed to calculate a derivative")]
	FiniteDifferenceError {
		#[from] source: FiniteDifferenceError
	}
}

// TODO: Check name convention
fn singlesteepdescent<T: Float + Debug + FromPrimitive + AddAssign, const LENGTH: usize>(func: fn([T; LENGTH]) -> T, x0: [T; LENGTH], lambda: T, h: [T; LENGTH], n: u8) -> Result<[T; LENGTH], SteepDescentError> {
	let dx: [T; LENGTH] = multivariate_central_finite_difference(func, x0, h, n)?;

	// TODO: Determine lambda (see document
	// TODO: Replace this with a Matrix type or variant
	let mut x1: [T; LENGTH] = [T::from(0).unwrap(); LENGTH];
	for i in 0..LENGTH {
		x1[i] = x0[i] - lambda * dx[i];
	}

	Ok(x1)
}

// TODO: Check max iter type and order of inputs to follow convention
pub fn steepdescent<T: Float + Debug + FromPrimitive + AddAssign, const LENGTH: usize>(func: fn([T; LENGTH]) -> T, x0: [T; LENGTH], lambda: T, h: [T; LENGTH], n: u8, max_iter: u32, toler: T) -> Result<[T; LENGTH], SteepDescentError> {
	let mut xn = x0;

	for _ in 0..max_iter {
		xn = singlesteepdescent(func, x0, lambda, h, n)?

		// TODO: How to measure tolerance with an array of values?
	}

	Ok(xn)
}

#[cfg(test)]
mod test {
	use num_traits::Pow;
	use super::steepdescent;

	fn case1(x: [f64; 2]) -> f64 {
		x[0] - x[1] + 2.0 * x[0].powi(2) + 2.0 * x[0] * x[1] + x[0].powi(2)
	}

	#[test]
	fn test_steepdescent() {

		assert_eq!(steepdescent(case1, [0.0, 0.0], 1.0, [1e-3; 2], 1, 1000, 1e-5).unwrap(), [-1.0, 1.5]);
	}
}