#[macro_use]
mod shared_ops;

#[cfg(backend = "64")]
pub(in crate::float_ops) mod f64_ops;
#[cfg(backend = "64")]
pub(in crate::float_ops) use f64_ops as backend;

#[cfg(backend = "rug")]
pub(in crate::float_ops) mod rug_ops;
#[cfg(backend = "rug")]
pub(in crate::float_ops) use rug_ops as backend;

#[cfg(backend = "32")]
pub(in crate::float_ops) mod f32_ops;
#[cfg(backend = "32")]
pub(in crate::float_ops) use f32_ops as backend;

// Special function algorithms shared by f32/f64 via SpecFloat trait.
// Rug uses native MPFR implementations instead.
#[cfg(any(backend = "64", backend = "32"))]
pub(in crate::float_ops) mod special;
