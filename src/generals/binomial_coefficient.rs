use std::fmt::Debug;
use num_traits::PrimInt;
use thiserror::Error;

#[derive(Error, Debug, Copy, Clone)]
pub enum FactorialError {
	#[error("The input to a factorial must be positive or 0")]
	InputMustBePositive
}


pub fn factorial<T: PrimInt + Debug>(n: T) -> Result<T, FactorialError>{
	if n < T::from(0).unwrap() {
		return Err(FactorialError::InputMustBePositive)
	}

	if n == T::from(0).unwrap() || n == T::from(1).unwrap() {
		return Ok(T::from(1.0).unwrap())
	}

	Ok(n * factorial(n - T::from(1.0).unwrap())?)
}

#[derive(Error, Debug, Copy, Clone)]
pub enum BinomialCoefficientError {
	#[error("The value of N must be larger than the value of K")]
	NmustBeLargest,
	#[error("One of the values is negative and therefore the factorial can not be calculated")]
	FactorialError{
		#[from] source: FactorialError
	},
}

pub fn binomial_coeff<T: PrimInt + Debug>(n: T, k: T) -> Result<T, BinomialCoefficientError> {
	if k > n {
		return Err(BinomialCoefficientError::NmustBeLargest)
	}

	Ok(factorial(n)? / ( factorial(k)? * factorial(n - k)?))
}