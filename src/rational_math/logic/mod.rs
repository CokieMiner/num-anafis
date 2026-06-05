#[cfg(any(backend = "64", backend = "32"))]
pub(in crate::rational_math) mod primitive_math;

#[cfg(any(backend = "64", backend = "32"))]
pub(in crate::rational_math) use primitive_math as backend;

// Override: rug (GMP-based arbitrary precision)
#[cfg(backend = "rug")]
pub(in crate::rational_math) mod rug_ops;

#[cfg(backend = "rug")]
pub(in crate::rational_math) use rug_ops as backend;
