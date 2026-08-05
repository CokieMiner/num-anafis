//! Numeric types with automatic representation promotion (Int → Rational → Float),
//! geometric algebra (Clifford) multivectors, and pluggable numeric backends (`backend32`,
//! `backend64`, `backendrug`).
//!
//! # Quick start
//! ```
//! use num_anafis::{Scalar, s, r, Number};
//! let a = s(3);                    // Scalar::Int(3)
//! let b = r(1, 3);                 // Scalar::Rational(1/3)
//! let c = a + b;                   // promoted to Rational: Scalar::Rational(10/3)
//! println!("{c}");                 // 10/3
//! ```
//!
//! Enable Clifford geometric algebra with the `clifford` feature:
//! ```
//! #[cfg(feature = "clifford")]
//! {
//! use num_anafis::{CliffordNumber, GeneratorSet, s};
//! let mv = CliffordNumber::scalar(GeneratorSet::empty(), s(2)).unwrap();
//! println!("{mv}");
//! }
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

// Backend features are mutually exclusive — enforce at compile time.
// (The build.rs falls back to backend64 if none is selected, so only
// the "multiple backends" case needs guarding here.)
#[cfg(any(
    all(feature = "backend32", feature = "backend64"),
    all(feature = "backend32", feature = "backendrug"),
    all(feature = "backend64", feature = "backendrug"),
))]
compile_error!(
    "Multiple `num-anafis` backends enabled simultaneously. \
     This usually happens when two dependencies hardcode different backends. \
     Library crates should use `default-features = false` and forward backend \
     features — see https://github.com/CokieMiner/num-anafis#backend-selection"
);

pub use error::NumAnafisError;
pub use float_ops::FloatType;
pub use int_math::IntType;
pub use number::{AnafisMathExt, IntoScalar, r, s};
pub use rational_math::RationalType;
pub use scalar::Scalar;
pub use traits::Number;

#[cfg(feature = "clifford")]
pub use clifford::{
    CliffordNumber, FastClifford, GeneratorSet, cga_gens, ci, e_minus, e_plus, e1, e2, e3, eps,
    inf, orig, pseudo3d, pseudo5d, qi, qj, qk, sj,
};

mod error;
mod float_ops;
mod int_math;
mod number;
mod rational_math;
mod scalar;
mod traits;
/// Complex numbers module.
pub mod complex;

#[cfg(feature = "clifford")]
mod clifford;

#[cfg(feature = "python")]
#[allow(
    missing_docs,
    reason = "Python bindings are not the main focus of this crate, so we can afford to be a bit less strict here."
)]
pub mod python_binding;
