#[macro_use]
mod shared_ops;

#[cfg(backend = "64")]
pub(in crate::int_math) mod i64_math;
#[cfg(backend = "64")]
pub(in crate::int_math) use i64_math as backend;

#[cfg(backend = "rug")]
pub(in crate::int_math) mod rug_int;
#[cfg(backend = "rug")]
pub(in crate::int_math) use rug_int as backend;

#[cfg(backend = "32")]
pub(in crate::int_math) mod i32_math;
#[cfg(backend = "32")]
pub(in crate::int_math) use i32_math as backend;
