#[macro_use]
mod shared_ops;

#[cfg(feature = "backend64")]
pub(in crate::int_math) mod i64_math;
#[cfg(feature = "backend64")]
pub(in crate::int_math) use i64_math as backend;

#[cfg(feature = "backendrug")]
pub(in crate::int_math) mod rug_int;
#[cfg(feature = "backendrug")]
pub(in crate::int_math) use rug_int as backend;

#[cfg(feature = "backend32")]
pub(in crate::int_math) mod i32_math;
#[cfg(feature = "backend32")]
pub(in crate::int_math) use i32_math as backend;
