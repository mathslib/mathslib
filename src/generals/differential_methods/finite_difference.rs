use std::fmt::Debug;
use std::ops::{Add, AddAssign};
use num_traits::{Float, Pow};
use crate::generals::{binomial_coeff, BinomialCoefficientError};

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


pub fn forward_finite_difference<T: Float + Debug + AddAssign>(func: fn(T) -> T, x0:T,  h: T, n: u8) -> Result<T, FiniteDifferenceError> {

	if h < T::from(0).unwrap(){
		return Err(FiniteDifferenceError::InvalidH)
	}

	let mut total: T = T::from(0).unwrap();

	for i in 0..(n+1) {
		total += T::from((-1.0).powi((n as i32 - i as i32))).unwrap() * T::from(binomial_coeff(n, i)?).unwrap() * func(x0 + T::from(i).unwrap() * h)
	}

	Ok(total / h)
}

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
