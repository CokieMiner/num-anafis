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

#![cfg_attr(
    all(
        not(feature = "std"),
        any(feature = "backend32", feature = "backend64", feature = "backendrug")
    ),
    no_std
)]

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
extern crate alloc;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
pub use error::NumAnafisError;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
pub use number::{AnafisMathExt, IntoScalar, r, s};

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
pub use traits::Number;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
pub use scalar::Scalar;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
pub use int_math::IntType;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
pub use float_ops::FloatType;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
pub use rational_math::RationalType;

#[cfg(feature = "clifford")]
pub use clifford::{
    CliffordNumber, FastClifford, GeneratorSet, cga_gens, ci, e_minus, e_plus, e1, e2, e3, eps,
    inf, orig, pseudo3d, pseudo5d, qi, qj, qk, sj,
};

// Backend features are mutually exclusive — enforce that exactly one is selected at compile time.
#[cfg(not(any(feature = "backend32", feature = "backend64", feature = "backendrug")))]
compile_error!(
    "Exactly one backend feature must be enabled: choose 'backend32', 'backend64', or 'backendrug'."
);

#[cfg(any(
    all(feature = "backend32", feature = "backend64"),
    all(feature = "backend32", feature = "backendrug"),
    all(feature = "backend64", feature = "backendrug"),
))]
compile_error!(
    "Multiple backend features cannot be enabled simultaneously: choose exactly one of 'backend32', 'backend64', or 'backendrug'."
);

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
mod error;

/// Core number types and operations — umbrella module for all numeric backends.
#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
mod number;

/// Shared traits used across numeric types.
#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
mod traits;

/// Generic scalar type with automatic representation management (Int → Rational → Float).
#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
mod scalar;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
mod int_math;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
mod float_ops;

#[cfg(any(feature = "backend32", feature = "backend64", feature = "backendrug"))]
mod rational_math;

#[cfg(feature = "clifford")]
mod clifford;

#[cfg(feature = "python")]
#[allow(
    missing_docs,
    reason = "Python bindings are not the main focus of this crate, so we can afford to be a bit less strict here."
)]
pub mod python_binding;
