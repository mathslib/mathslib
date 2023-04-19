use std::fmt::Debug;
use std::ops::AddAssign;
use num_traits::Float;
use crate::generals::{binomial_coeff, BinomialCoefficientError};
use thiserror::Error;

// TODO: Force h to be greater than 0
// TODO: Consider grouping the similar logic into the same function

#[derive(Error, Debug)]
pub enum FiniteDifferenceError {
	#[error("Invalid value of h")]
	InvalidH,
	#[error("Failed to calculate the binomial coefficient for the order of the FDM.")]
	BinomialCoefficientError {
		#[from] source: BinomialCoefficientError
	}
}

/// Finite difference methods are the simplest, more well known methods for numerically determining an approximation for the derivative of a function.
///
/// The [central finite difference](central_finite_difference) tends to be a more accurate approximation.
///
/// # Arguments
/// * func: The function to determine the derivative for
/// * x0: The x value to calculate the derivative for
/// * h: The spacing of the bound considered (A smaller value will give a more accurate result but caution must be taken to not loose resolution)
/// * n: Order of the finite difference method. Higher value get better approximations but with a larger number of evaluations required.
///
/// # Example
/// ```
/// use mathslib::generals::differential_methods::forward_finite_difference;
///
/// fn case1(x: f64) -> f64 {
/// 	x.powi(2) + 6.0 * x + 3.0
/// }
///
/// fn main(){
/// 	assert_eq!(forward_finite_difference(case1, 3.0, 1e-5, 1).unwrap().round(), 12.0)
/// }
/// ```
///
pub fn forward_finite_difference<T: Float + Debug + AddAssign>(func: fn(T) -> T, x0:T,  h: T, n: u8) -> Result<T, FiniteDifferenceError> {

	if h < T::from(0).unwrap(){
		return Err(FiniteDifferenceError::InvalidH)
	}

	let mut total: T = T::from(0).unwrap();

	for i in 0..(n+1) {
		total += T::from((-1.0).powi(n as i32 - i as i32)).unwrap() * T::from(binomial_coeff(n, i)?).unwrap() * func(x0 + T::from(i).unwrap() * h)
	}

	Ok(total / h)
}

/// Finite difference methods are the simplest, more well known methods for numerically determining an approximation for the derivative of a function.
///
/// The [central finite difference](central_finite_difference) tends to be a more accurate approximation.
///
/// # Arguments
/// * func: The function to determine the derivative for
/// * x0: The x value to calculate the derivative for
/// * h: The spacing of the bound considered (A smaller value will give a more accurate result but caution must be taken to not loose resolution)
/// * n: Order of the finite difference method. Higher value get better approximations but with a larger number of evaluations required.
///
/// # Example
/// ```
/// use mathslib::generals::differential_methods::backwards_finite_difference;
///
/// fn case1(x: f64) -> f64 {
/// 	x.powi(2) + 6.0 * x + 3.0
/// }
///
/// fn main(){
/// 	assert_eq!(backwards_finite_difference(case1, 3.0, 1e-5, 1).unwrap().round(), 12.0)
/// }
/// ```
///
pub fn backwards_finite_difference<T: Float + Debug + AddAssign>(func: fn(T) -> T, x0: T, h: T, n: u8) -> Result<T, FiniteDifferenceError> {

	if h < T::from(0).unwrap(){
		return Err(FiniteDifferenceError::InvalidH)
	}

	let mut total: T = T::from(0).unwrap();

	for i in 0..(n+1) {
		total += T::from((-1.0).powi(i as i32)).unwrap() * T::from(binomial_coeff(n, i)?).unwrap() * func(x0 - T::from(i).unwrap() * h)
	}

	Ok(total/h)
}

/// Finite difference methods are the simplest, more well known methods for numerically determining an approximation for the derivative of a function.
///
/// This is the most accurate finite difference method with an estimated truncation error in the order of O(h^2)
///
/// # Arguments
/// * func: The function to determine the derivative for
/// * x0: The x value to calculate the derivative for
/// * h: The spacing of the bound considered (A smaller value will give a more accurate result but caution must be taken to not loose resolution)
/// * n: Order of the finite difference method. Higher value get better approximations but with a larger number of evaluations required.
///
/// # Example
/// ```
/// use mathslib::generals::differential_methods::central_finite_difference;
///
/// fn case1(x: f64) -> f64 {
/// 	x.powi(2) + 6.0 * x + 3.0
/// }
///
/// fn main(){
/// 	assert_eq!(central_finite_difference(case1, 3.0, 1e-5, 1).unwrap().round(), 12.0)
/// }
/// ```
///
pub fn central_finite_difference<T: Float + Debug + AddAssign>(func: fn(T) -> T, x0: T, h: T, n: u8) -> Result<T, FiniteDifferenceError> {

	if h < T::from(0).unwrap(){
		return Err(FiniteDifferenceError::InvalidH)
	}

	let mut total: T = T::from(0).unwrap();

	for i in 0..(n+1){
		total += T::from((-1.0).powi(i as i32)).unwrap() * T::from(binomial_coeff(n, i)?).unwrap() * func(x0 + ((T::from(n).unwrap() / T::from(2).unwrap()) - T::from(i).unwrap()) * h)
	}

	Ok(total/h)
}
