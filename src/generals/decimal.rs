use std::fmt::Debug;
use num_traits::{Float, FromPrimitive};

pub trait Decimal {

	/// Round a decimal number to the specified number of decimal places
	///
	/// # Arguments
	/// * dp: Number of decimal places
	fn round_dp(&self, dp: i32) -> Self;
}

impl<T: Float + PartialOrd + Debug + FromPrimitive> Decimal for T {
	fn round_dp(&self, dp: i32) -> T{
		let ten = T::from_f64(10.0).unwrap().powi(dp);
		let a = *self * ten;
		a.round() / ten
	}
}
