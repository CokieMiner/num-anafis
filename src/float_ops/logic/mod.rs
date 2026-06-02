#[macro_use]
mod shared_ops;

#[cfg(feature = "backend64")]
pub(in crate::float_ops) mod f64_ops;
#[cfg(feature = "backend64")]
pub(in crate::float_ops) use f64_ops as backend;

#[cfg(feature = "backendrug")]
pub(in crate::float_ops) mod rug_ops;
#[cfg(feature = "backendrug")]
pub(in crate::float_ops) use rug_ops as backend;

#[cfg(feature = "backend32")]
pub(in crate::float_ops) mod f32_ops;
#[cfg(feature = "backend32")]
pub(in crate::float_ops) use f32_ops as backend;

// Special function algorithms shared by f32/f64 via SpecFloat trait.
// Rug uses native MPFR implementations instead.
#[cfg(any(feature = "backend64", feature = "backend32"))]
pub(in crate::float_ops) mod special;
