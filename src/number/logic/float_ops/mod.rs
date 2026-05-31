#[cfg(feature = "backend64")]
mod f64_ops;

#[cfg(feature = "backendrug")]
mod rug_ops;

#[cfg(feature = "backend32")]
mod f32_ops;

// Special function algorithms shared by f32/f64 via SpecFloat trait.
// Rug uses native MPFR implementations instead.
#[cfg(any(feature = "backend64", feature = "backend32"))]
pub(super) mod special;

mod api;

pub use api::*;
