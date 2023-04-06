#[cfg(test)]
mod scalar_optimization {
	use num_traits::Pow;
	use crate::generals::Decimal;
	use crate::optimize::scalar::bracket_optimizers::{bracket_gr_minimize};
	use crate::optimize::scalar::bound_optimizers::{bound_gr_minimize as bound_gr};

	fn case_1(x: f64) -> f64{x.pow(2) + 6.0*x + 3.0	}

	#[test]
	fn golden_ratio_case_1() {
		assert_eq!(bracket_gr_minimize::<f64>(case_1, 4.0, -9.0, 1.0, 1e-4, 2000).unwrap().round_dp(4), -3.0000)
	}

	#[test]
	fn golden_ratio_case_1_bound_method() {
		assert_eq!(bound_gr::<f64>(case_1, 4.0, -9.0, 1e-4, 2000).unwrap().round_dp(4), -3.0000)
	}

	fn case_2(x: f64) -> f64 {x.pow(4) - 2.0 * x.pow(3) + 4.0}

	#[test]
	fn golden_ratio_case_2() {
		assert_eq!(bracket_gr_minimize::<f64>(case_2, 1.0, 2.0, 0.5, 1e-4, 2000).unwrap().round_dp(4), 1.5000)
	}

	#[test]
	fn golden_ratio_case_2_bound_method() {
		assert_eq!(bound_gr::<f64>(case_2, 0.5, 2.0, 1e-4, 2000).unwrap().round_dp(4), 1.5000)
	}
}