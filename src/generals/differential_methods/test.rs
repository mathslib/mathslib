use crate::generals::{binomial_coeff, factorial};
use super::forward_finite_difference;

fn case1(x: f64) -> f64 {
	x.powi(2) + 6.0 * x + 3.0
}

#[test]
fn forward_finite_difference_test(){
	assert_eq!(forward_finite_difference(case1, 3.0, 1e-5, 1).unwrap().round(), 12.0)
}
