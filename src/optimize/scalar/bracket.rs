use num_traits::{abs, Float, FromPrimitive, Pow, Signed};
use std::cmp::Ordering::Equal;
use std::fmt::Debug;
use thiserror::Error;

#[derive(Debug, Copy, Clone)]
struct Bracket<T: Float + PartialOrd + Debug> {
    pub left: T,
    pub f_left: T,
    pub right: T,
    pub f_right: T,
    pub center: T,
    pub f_center: T,
}

#[derive(Error, Debug, Copy, Clone)]
pub enum BracketError {
    #[error("A bracket bound is duplicated and therefore a bracket would be invalid")]
    DupeBoundForBracket,
    #[error("The specified bracket does not encapsulate a dip. It might have a local minima within the bound but this is unknown")]
    BracketNotADip,
}

impl<T: Float + PartialOrd + Debug> Bracket<T> {
    pub fn new(x1: T, x2: T, x3: T, func: fn(T) -> T) -> Result<Self, BracketError> {
        let mut elements: [T; 3] = [x1, x2, x3];

        elements.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Equal));

        Self::new_processed(
            elements[0],
            func(elements[0]),
            elements[1],
            func(elements[1]),
            elements[2],
            func(elements[2]),
        )
    }

    pub fn new_processed(
        left: T,
        f_left: T,
        center: T,
        f_center: T,
        right: T,
        f_right: T,
    ) -> Result<Self, BracketError> {
        if left == center || center == right {
            return Err(BracketError::DupeBoundForBracket);
        }

        if !((f_center < f_left) && (f_center < f_right)) {
            return Err(BracketError::BracketNotADip);
        }

        Ok(Bracket {
            left,
            f_left,
            right,
            f_right,
            center,
            f_center,
        })
    }

    pub fn longer_bound(&self, func: fn(T) -> T, new_val: T) -> Result<Bracket<T>, BracketError> {
        let f_new_val = func(new_val);

        let new: Self;

        if (f_new_val < self.f_center) & (new_val < self.center) {
            new = Self::new_processed(
                self.left,
                self.f_left,
                new_val,
                f_new_val,
                self.center,
                self.f_center,
            )?;
        } else if (f_new_val > self.f_center) & (new_val < self.center) {
            new = Self::new_processed(
                new_val,
                f_new_val,
                self.center,
                self.f_center,
                self.right,
                self.f_right,
            )?;
        } else if (f_new_val < self.f_center) & (new_val > self.center) {
            new = Self::new_processed(
                self.center,
                self.f_center,
                new_val,
                f_new_val,
                self.right,
                self.f_right,
            )?;
        } else {
            new = Self::new_processed(
                self.left,
                self.f_left,
                self.center,
                self.f_center,
                new_val,
                f_new_val,
            )?;
        }

        Ok(new)
    }

    pub fn new_val_from_ratio(&self, ratio: T) -> T {
        if self.right - self.center > self.center - self.left {
            return self.right - (self.right - self.center) / ratio;
        }
        self.left + (self.center - self.left) / ratio
    }

    pub fn parabolic_interpolation(&self) -> T {
        ((self.f_right - self.f_center) * (self.center.powi(2) - self.left.powi(2))
            + (self.f_left - self.f_center) * (self.right.powi(2) - self.center.powi(2)))
            / (T::from(2.0).unwrap()
                * ((self.f_right - self.f_center) * (self.center - self.left)
                    + (self.f_left - self.f_center) * (self.right - self.center)))
    }
}

#[derive(Error, Debug, Copy, Clone)]
pub enum BracketRatioOptimizerError {
    #[error("A bracket is invalid")]
    BracketError {
        #[from]
        source: BracketError,
    },
    #[error("Invalid ratio, the ratio must be greater than 1 and smaller than 2")]
    InvalidRatio,
    #[error("Invalid tolerance must be greater than 0")]
    InvalidTolerance,
}

#[allow(dead_code)]
/// The iterative bracket based approach for minimization of a problem using the golden ratio.
///
/// It is recommended to use the [bound method](bound_gr_minimize) as it tends to converge faster than this method.
///
/// # Arguments
/// * func: The objective function to minimize. The function, of the form f(x: T) -> T, must take in a single scalar value and return another.
/// * xi: The three values characterizing the bracket
/// * tolerance: The tolerance requirement to determine convergence
/// * max_iter: The maximum number of iterations to loop over.
///
/// # Example
/// ```
/// use mathslib::optimize::scalar::bracket_optimizers::bracket_gr_minimize;
/// use mathslib::generals::Decimal;
///
/// fn case_1(x: f64) -> f64{x*x + 6.0*x + 3.0	}
///
/// fn main() {
/// 	assert_eq!(bracket_gr_minimize::<f64>(case_1, 4.0, -9.0, 1.0, 1e-4, 2000).unwrap().round_dp(4), -3.0000)
///}
/// ```
pub fn bracket_gr_minimize<T: Float + FromPrimitive + Signed + PartialOrd + Debug>(
    func: fn(T) -> T,
    x1: T,
    x2: T,
    x3: T,
    tolerance: T,
    max_iter: u32,
) -> Result<T, BracketRatioOptimizerError> {
    bracket_ratio_minimize(
        func,
        x1,
        x2,
        x3,
        T::from_f64((5.0.pow(0.5) + 1.0) / 2.0).unwrap(),
        tolerance,
        max_iter,
    )
}

fn single_bracket_minimize<T: Float + PartialOrd + Debug + FromPrimitive>(
    func: fn(T) -> T,
    bounds: Bracket<T>,
    ratio: T,
) -> Result<Bracket<T>, BracketRatioOptimizerError> {
    // Note this function will never return an invalid tolerance error. Typically, a different error struct would be generated but, because it is a private function, it should never impact the library user.

    if ratio < T::from_f64(1.0).unwrap() {
        return Err(BracketRatioOptimizerError::InvalidRatio);
    }

    Ok(bounds.longer_bound(func, bounds.new_val_from_ratio(ratio))?)
}

#[allow(dead_code)]
/// The iterative bracket based approach for minimization of a problem.
///
/// It is recommended to use the [bound method](bound_minimize) as it tends to converge faster and is less error prone than this method.
///
/// # Arguments
/// * func: The objective function to minimize. The function, of the form f(x: T) -> T, must take in a single scalar value and return another.
/// * xi: The three values characterizing the bracket
/// * ratio: The ratio around a bound to split around
/// * tolerance: The tolerance requirement to determine convergence
/// * max_iter: The maximum number of iterations to loop over.
///
/// # Example
/// ```
/// use mathslib::optimize::scalar::bracket_optimizers::bracket_ratio_minimize;
/// use mathslib::generals::Decimal;
///
/// fn case_1(x: f64) -> f64{x*x + 6.0*x + 3.0	}
///
/// fn main() {
/// 	assert_eq!(bracket_ratio_minimize::<f64>(case_1, 4.0, -9.0, 1.0, 1.5, 1e-4, 2000).unwrap().round_dp(4), -3.0000)
///}
/// ```
pub fn bracket_ratio_minimize<T: Float + FromPrimitive + Signed + PartialOrd + Debug>(
    func: fn(T) -> T,
    x1: T,
    x2: T,
    x3: T,
    ratio: T,
    tolerance: T,
    max_iter: u32,
) -> Result<T, BracketRatioOptimizerError> {
    if tolerance < T::from_f64(0.0).unwrap() {
        return Err(BracketRatioOptimizerError::InvalidTolerance);
    }

    let mut bounds = Bracket::new(x1, x2, x3, func)?;
    let mut old_val: T = bounds.f_center + tolerance;

    for _ in 0..max_iter {
        bounds = match &single_bracket_minimize(func, bounds, ratio) {
            Ok(bracket) => *bracket,
            Err(error) => {
                match error {
                    // A duplicated value in the bracket typically means that we have reached the maximum level of precision posible by the previous iteration
                    // FIXME: Consider if this should return the error or assume it has converge to the best possible value it could
                    BracketRatioOptimizerError::BracketError { .. } => break,
                    _ => return Err(*error),
                }
            }
        };

        if abs(old_val - bounds.f_center) < tolerance {
            break;
        }

        old_val = bounds.center;
    }

    Ok(bounds.center)
}

#[allow(dead_code)]
/// Perform parabolic interpolation over a bracket to find an approximation to the minimum.
///
/// # Example
///
/// ```
/// use mathslib::optimize::scalar::bracket_optimizers::parabolic_interpolation;
///
/// fn case(x: f64) -> f64{
///     x.powi(4) - 2.0 * x.powi(3) + 4.0
/// }
///
/// assert_eq!(parabolic_interpolation(case, 0.5, 1.0, 2.0).unwrap(), 1.2142857142857142)
///
/// ```
///
pub fn parabolic_interpolation<T: Float + Debug>(
    func: fn(T) -> T,
    x1: T,
    x2: T,
    x3: T,
) -> Result<T, BracketError> {
    Ok(Bracket::new(x1, x2, x3, func)?.parabolic_interpolation())
}

#[derive(Error, Debug, Copy, Clone)]
pub enum BracketPIOptimizerError {
    #[error("A bracket is invalid")]
    BracketError {
        #[from]
        source: BracketError,
    },
    #[error("Invalid tolerance, must be greater than 0")]
    InvalidTolerance,
    #[error("Reached a dead end")]
    DeadEnd,
}

/// Perform parabolic interpolation until converging to a minima.
///
/// # Arguments
/// * func: The objective function to minimize. The function, of the form f(x: T) -> T, must take in a single scalar value and return another.
/// * xi: The values forming the bracket
/// * tolerance: The tolerance requirement to determine convergence
/// * max_iter: The maximum number of iterations to loop over.
///
/// # Example
/// ```
/// 	use mathslib::optimize::scalar::bracket_optimizers::bracket_pi_minimize;
/// 	use mathslib::generals::Decimal;
/// 	use num_traits::Pow;
///
///     fn case_2(x: f64) -> f64 {
/// 		x.pow(4) - 2.0 * x.pow(3) + 4.0
/// 	}
///
///     fn main() {
///
/// 		assert_eq!(bracket_pi_minimize(case_2, 0.5, 1.0, 2.0, 1e-4, 5).unwrap().round_dp(6), 1.482046)
///     }
/// ```
///
pub fn bracket_pi_minimize<T: Float + Debug + FromPrimitive + Signed>(
    func: fn(T) -> T,
    x1: T,
    x2: T,
    x3: T,
    tolerance: T,
    max_iter: u32,
) -> Result<T, BracketPIOptimizerError> {
    // Validate the tolerance
    if tolerance < T::from_f64(0.0).unwrap() {
        return Err(BracketPIOptimizerError::InvalidTolerance);
    }

    // Generate the initial bracket
    let mut bounds = Bracket::new(x1, x2, x3, func)?;
    let mut old_val: T = bounds.f_center + tolerance;

    for _ in 0..max_iter {
        bounds = bounds.longer_bound(func, bounds.parabolic_interpolation())?;

        /// Check if reached a dead end (Dead ends is when a bracket results in the same bracket
        /// An example would be, given a bracket [(0, 6), (1, 4), (2, 6)] which results in the same bracket but not necessarily a minimum
        if old_val == bounds.f_center {
            return Err(BracketPIOptimizerError::DeadEnd);
        }

        if abs(old_val - bounds.f_center) < tolerance {
            break;
        }

        old_val = bounds.center
    }

    Ok(bounds.center)
}
