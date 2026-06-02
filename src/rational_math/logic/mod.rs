#[cfg(any(feature = "backend64", feature = "backend32"))]
pub(in crate::rational_math) mod primitive_math;

#[cfg(any(feature = "backend64", feature = "backend32"))]
pub(in crate::rational_math) use primitive_math as backend;

// Override: rug (GMP-based arbitrary precision)
#[cfg(feature = "backendrug")]
pub(in crate::rational_math) mod rug_ops;

#[cfg(feature = "backendrug")]
pub(in crate::rational_math) use rug_ops as backend;
