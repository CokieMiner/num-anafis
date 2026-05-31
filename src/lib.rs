//! Numeric core for `SymbAnaFis`.

#![cfg_attr(
    all(
        not(feature = "std"),
        any(feature = "backend32", feature = "backend64", feature = "backendrug")
    ),
    no_std
)]

#[cfg(all(any(feature = "backend32", feature = "backend64", feature = "backendrug")))]
extern crate alloc;

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

#[cfg(all(any(feature = "backend32", feature = "backend64", feature = "backendrug")))]
mod error;

#[cfg(all(any(feature = "backend32", feature = "backend64", feature = "backendrug")))]
mod number;

#[cfg(all(any(feature = "backend32", feature = "backend64", feature = "backendrug")))]
pub use error::NumAnafisError;

#[cfg(all(any(feature = "backend32", feature = "backend64", feature = "backendrug")))]
pub use number::{
    AnafisMathExt, FloatType, IntType, IntoScalar, Number, RationalType, Scalar, r, s,
};

#[cfg(all(feature = "clifford"))]
pub use number::{
    CliffordNumber, FastClifford, GeneratorSet, cga_gens, ci, e_minus, e_plus, e1, e2, e3, eps,
    inf, orig, pseudo3d, pseudo5d, qi, qj, qk, sj,
};

#[cfg(all(feature = "python"))]
#[allow(
    missing_docs,
    reason = "Python bindings are not the main focus of this crate, so we can afford to be a bit less strict here."
)]
pub mod python_binding;
