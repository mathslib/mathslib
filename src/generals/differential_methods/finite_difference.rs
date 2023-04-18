use std::fmt::Debug;
use std::ops::AddAssign;
use num_traits::{Float, Pow};
use crate::generals::{binomial_coeff, BinomialCoefficientError};


pub fn forward_finite_difference<T: Float + Debug + AddAssign>(func: fn(T) -> T, x0:T,  h: T, n: u8) -> Result<T, BinomialCoefficientError> {

	let mut total: T = T::from(0).unwrap();

	for i in 0..(n+1) {
		total += T::from((-1.0).powi((n as i32 - i as i32))).unwrap() * T::from(binomial_coeff(n, i)?).unwrap() * func(x0 + T::from(i).unwrap() * h)
	}

	return Ok(total / h)
}

