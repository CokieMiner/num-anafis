use alloc::vec::Vec;
use core::array;
use core::cmp::Ordering;
use core::fmt::{self, Display, Formatter};
use core::hash::{Hash, Hasher};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::types::{CliffordNumber, GeneratorSet, INLINE_COEFF_COUNT, scalar_zero};
use crate::scalar::Scalar;
use crate::traits::Number;
use crate::error::NumAnafisError;

/// A purely algebraic, zero-allocation Clifford multivector optimized for `#![no_std]`
/// and const-generic compile-time evaluation.
///
/// `P`, `Q`, `R` define the signature (positive, negative, zero metrics).
/// Currently limits to algebras where `P + Q + R <= 5` for inline storage.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(::serde::Serialize, ::serde::Deserialize))]
#[non_exhaustive]
pub struct FastClifford<const P: usize, const Q: usize, const R: usize> {
    /// The array of multivector coefficients, strictly bounded to `INLINE_COEFF_COUNT`.
    pub coeffs: [Scalar; INLINE_COEFF_COUNT],
}

impl<const P: usize, const Q: usize, const R: usize> FastClifford<P, Q, R> {
    const _CHECK_SIZE: () = assert!(
        P + Q + R <= 5,
        "FastClifford requires P + Q + R <= 5 for inline storage"
    );

    /// Creates a multivector initialized completely to zero.
    #[must_use]
    pub fn zero() -> Self {
        let () = Self::_CHECK_SIZE;
        Self {
            coeffs: array::from_fn(|_| scalar_zero()),
        }
    }
    /// Sets to zero any coefficient whose absolute value is strictly less than `Scalar::epsilon()`.
    ///
    /// This is useful to clean up `$10^{-16}$` numerical noise.
    #[must_use]
    pub fn chop(mut self) -> Self {
        let eps = Scalar::epsilon();
        for c in &mut self.coeffs {
            if c.abs().total_cmp(&eps) == Ordering::Less {
                *c = scalar_zero();
            }
        }
        self
    }
    /// Returns the total number of active generators mapped to this algebra (`P + Q + R`).
    #[must_use]
    pub const fn active_generators() -> usize {
        P + Q + R
    }
    /// Evaluates the geometric metric for generator index `i` (+1, -1, or 0).
    ///
    /// # Panics
    /// Panics if the requested index `i` is greater than or equal to `P + Q + R`.
    #[must_use]
    #[allow(
        clippy::panic,
        reason = "Out-of-bounds metric access is a fundamental logic error, equivalent to an array out-of-bounds index. It is unrecoverable mathematically."
    )]
    pub const fn metric_at(i: usize) -> i8 {
        if i < P {
            1
        } else if i < P + Q {
            -1
        } else if i < P + Q + R {
            0
        } else {
            panic!("Generator index out of bounds")
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> From<&FastClifford<P, Q, R>>
    for CliffordNumber
{
    #[allow(
        clippy::cast_possible_truncation,
        reason = "i < P+Q+R <= 5, always fits in u32"
    )]
    fn from(val: &FastClifford<P, Q, R>) -> Self {
        let n = P + Q + R;
        let mut entries = Vec::with_capacity(n);
        for i in 0..P {
            entries.push((i as u32, 1));
        }
        for i in P..P + Q {
            entries.push((i as u32, -1));
        }
        for i in P + Q..n {
            entries.push((i as u32, 0));
        }
        let gens = GeneratorSet::from_sorted(entries);

        let mut mv = Self::zero_unchecked(gens);
        let len = 1 << n;
        for (i, coeff) in val.coeffs.iter().enumerate().take(len) {
            mv.coeffs_mut_slice()[i] = coeff.clone();
        }
        mv
    }
}

impl<const P: usize, const Q: usize, const R: usize> FastClifford<P, Q, R> {
    /// Creates a `FastClifford` from a `CliffordNumber` without checking generator match.
    #[must_use]
    pub(crate) fn from_unchecked(val: &CliffordNumber) -> Self {
        let () = Self::_CHECK_SIZE;
        let mut coeffs = array::from_fn(|_| scalar_zero());
        let len = 1 << (P + Q + R);
        for (i, coeff) in coeffs.iter_mut().enumerate().take(len) {
            *coeff = val.coeff(i).clone();
        }
        Self { coeffs }
    }
}

impl<const P: usize, const Q: usize, const R: usize> TryFrom<&CliffordNumber>
    for FastClifford<P, Q, R>
{
    type Error = NumAnafisError;

    fn try_from(val: &CliffordNumber) -> Result<Self, Self::Error> {
        let () = Self::_CHECK_SIZE;
        let gens = val.generator_set();
        if gens.len() != P + Q + R {
            return Err(NumAnafisError::MismatchedGeneratorSet);
        }
        for i in 0..P {
            if gens.metric_at(i) != 1 {
                return Err(NumAnafisError::MismatchedGeneratorSet);
            }
        }
        for i in P..P + Q {
            if gens.metric_at(i) != -1 {
                return Err(NumAnafisError::MismatchedGeneratorSet);
            }
        }
        for i in P + Q..P + Q + R {
            if gens.metric_at(i) != 0 {
                return Err(NumAnafisError::MismatchedGeneratorSet);
            }
        }
        Ok(Self::from_unchecked(val))
    }
}

impl<const P: usize, const Q: usize, const R: usize> Display for FastClifford<P, Q, R> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mv = CliffordNumber::from(self);
        write!(f, "{mv}")
    }
}

impl<const P: usize, const Q: usize, const R: usize> PartialOrd for FastClifford<P, Q, R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let mv1 = CliffordNumber::from(self);
        let mv2 = CliffordNumber::from(other);
        mv1.partial_cmp(&mv2)
    }
}

impl<const P: usize, const Q: usize, const R: usize> Hash for FastClifford<P, Q, R> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mv = CliffordNumber::from(self);
        mv.hash(state);
    }
}

macro_rules! impl_op {
    ($trait:ident, $method:ident) => {
        impl<const P: usize, const Q: usize, const R: usize> $trait for FastClifford<P, Q, R> {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self::Output {
                let mv1 = CliffordNumber::from(&self);
                let mv2 = CliffordNumber::from(&rhs);
                let res = mv1.$method(&mv2);
                Self::from_unchecked(&res)
            }
        }
        impl<'num, const P: usize, const Q: usize, const R: usize>
            $trait<&'num FastClifford<P, Q, R>> for FastClifford<P, Q, R>
        {
            type Output = Self;
            fn $method(self, rhs: &'num FastClifford<P, Q, R>) -> Self::Output {
                let mv1 = CliffordNumber::from(&self);
                let mv2 = CliffordNumber::from(rhs);
                let res = mv1.$method(&mv2);
                Self::from_unchecked(&res)
            }
        }
    };
}
// Direct zero-allocation implementations for simple coefficient-wise operations.
impl<const P: usize, const Q: usize, const R: usize> Add for FastClifford<P, Q, R> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let sum = &self.coeffs[i] + &rhs.coeffs[i];
            self.coeffs[i] = sum;
        }
        self
    }
}

impl<'num, const P: usize, const Q: usize, const R: usize> Add<&'num Self>
    for FastClifford<P, Q, R>
{
    type Output = Self;
    fn add(mut self, rhs: &'num Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let sum = &self.coeffs[i] + &rhs.coeffs[i];
            self.coeffs[i] = sum;
        }
        self
    }
}

impl<const P: usize, const Q: usize, const R: usize> Sub for FastClifford<P, Q, R> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let diff = &self.coeffs[i] - &rhs.coeffs[i];
            self.coeffs[i] = diff;
        }
        self
    }
}

impl<'num, const P: usize, const Q: usize, const R: usize> Sub<&'num Self>
    for FastClifford<P, Q, R>
{
    type Output = Self;
    fn sub(mut self, rhs: &'num Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let diff = &self.coeffs[i] - &rhs.coeffs[i];
            self.coeffs[i] = diff;
        }
        self
    }
}

// Geometric product and division still require the spectral/Cayley pipeline.
impl_op!(Mul, mul);
impl_op!(Div, div);

impl<const P: usize, const Q: usize, const R: usize> Neg for FastClifford<P, Q, R> {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let c = -&self.coeffs[i];
            self.coeffs[i] = c;
        }
        self
    }
}

// ============================================================================
// Compound assignment operators
// ============================================================================

impl<const P: usize, const Q: usize, const R: usize> AddAssign for FastClifford<P, Q, R> {
    fn add_assign(&mut self, rhs: Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let sum = &self.coeffs[i] + &rhs.coeffs[i];
            self.coeffs[i] = sum;
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> AddAssign<&Self> for FastClifford<P, Q, R> {
    fn add_assign(&mut self, rhs: &Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let sum = &self.coeffs[i] + &rhs.coeffs[i];
            self.coeffs[i] = sum;
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> SubAssign for FastClifford<P, Q, R> {
    fn sub_assign(&mut self, rhs: Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let diff = &self.coeffs[i] - &rhs.coeffs[i];
            self.coeffs[i] = diff;
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> SubAssign<&Self> for FastClifford<P, Q, R> {
    fn sub_assign(&mut self, rhs: &Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            let diff = &self.coeffs[i] - &rhs.coeffs[i];
            self.coeffs[i] = diff;
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> MulAssign for FastClifford<P, Q, R> {
    fn mul_assign(&mut self, rhs: Self) {
        let mv1 = CliffordNumber::from(&*self);
        let mv2 = CliffordNumber::from(&rhs);
        *self = Self::from_unchecked(&(mv1 * &mv2));
    }
}

impl<const P: usize, const Q: usize, const R: usize> MulAssign<&Self> for FastClifford<P, Q, R> {
    fn mul_assign(&mut self, rhs: &Self) {
        let mv1 = CliffordNumber::from(&*self);
        let mv2 = CliffordNumber::from(rhs);
        *self = Self::from_unchecked(&(mv1 * &mv2));
    }
}

impl<const P: usize, const Q: usize, const R: usize> DivAssign for FastClifford<P, Q, R> {
    fn div_assign(&mut self, rhs: Self) {
        let mv1 = CliffordNumber::from(&*self);
        let mv2 = CliffordNumber::from(&rhs);
        *self = Self::from_unchecked(&(mv1 / &mv2));
    }
}

impl<const P: usize, const Q: usize, const R: usize> DivAssign<&Self> for FastClifford<P, Q, R> {
    fn div_assign(&mut self, rhs: &Self) {
        let mv1 = CliffordNumber::from(&*self);
        let mv2 = CliffordNumber::from(rhs);
        *self = Self::from_unchecked(&(mv1 / &mv2));
    }
}

macro_rules! delegate_number_unary {
    ($($method:ident),*) => {
        $(
            fn $method(&self) -> Self {
                let mv = CliffordNumber::from(self);
                let res = mv.$method();
                Self::from_unchecked(&res)
            }
        )*
    };
}

macro_rules! delegate_number_binary {
    ($($method:ident),*) => {
        $(
            fn $method(&self, other: &Self) -> Self {
                let mv1 = CliffordNumber::from(self);
                let mv2 = CliffordNumber::from(other);
                let res = mv1.$method(&mv2);
                Self::from_unchecked(&res)
            }
        )*
    };
}

impl<const P: usize, const Q: usize, const R: usize> Number for FastClifford<P, Q, R> {
    delegate_number_unary!(
        sin, cos, tan, cot, sec, csc, asin, acos, atan, acot, asec, acsc, sinh, cosh, tanh, coth,
        sech, csch, asinh, acosh, atanh, acoth, acsch, asech, exp, expm1, exp_neg, ln, log1p, sqrt,
        cbrt, abs, signum, floor, ceil, round, fract, negate, erf, erfc, gamma, lgamma, digamma,
        trigamma, tetragamma, sinc, elliptic_k, elliptic_e, zeta, exp_polar
    );

    delegate_number_binary!(
        atan2, log_base, pow, besselj, bessely, besseli, besselk, polygamma, beta, zeta_deriv,
        lambertw, hermite, num_max, num_min
    );

    fn assoc_legendre(&self, l: &Self, m: &Self) -> Self {
        let mv1 = CliffordNumber::from(self);
        let mv_l = CliffordNumber::from(l);
        let mv_m = CliffordNumber::from(m);
        let res = mv1.assoc_legendre(&mv_l, &mv_m);
        Self::from_unchecked(&res)
    }

    fn spherical_harmonic(&self, l: &Self, m: &Self, phi: &Self) -> Self {
        let mv1 = CliffordNumber::from(self);
        let mv_l = CliffordNumber::from(l);
        let mv_m = CliffordNumber::from(m);
        let mv_phi = CliffordNumber::from(phi);
        let res = mv1.spherical_harmonic(&mv_l, &mv_m, &mv_phi);
        Self::from_unchecked(&res)
    }

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(Number::is_zero)
    }

    fn is_one(&self) -> bool {
        self.coeffs[0].is_one() && self.coeffs.iter().skip(1).all(Number::is_zero)
    }

    fn is_neg_one(&self) -> bool {
        self.coeffs[0].is_neg_one() && self.coeffs.iter().skip(1).all(Number::is_zero)
    }

    fn is_integer(&self) -> bool {
        let len = 1 << Self::active_generators();
        self.coeffs[0].is_integer()
            && self
                .coeffs
                .iter()
                .skip(1)
                .take(len - 1)
                .all(Number::is_zero)
    }

    fn is_negative(&self) -> bool {
        let len = 1 << Self::active_generators();
        self.coeffs
            .iter()
            .skip(1)
            .take(len - 1)
            .all(Number::is_zero)
            && self.coeffs[0].is_negative()
    }

    fn is_positive(&self) -> bool {
        let len = 1 << Self::active_generators();
        self.coeffs
            .iter()
            .skip(1)
            .take(len - 1)
            .all(Number::is_zero)
            && self.coeffs[0].is_positive()
    }

    fn is_finite(&self) -> bool {
        self.coeffs.iter().all(Number::is_finite)
    }

    fn to_float(&self) -> Self {
        let mut coeffs = array::from_fn(|_| scalar_zero());
        for (i, c) in self.coeffs.iter().enumerate().take(INLINE_COEFF_COUNT) {
            coeffs[i] = c.to_float();
        }
        Self { coeffs }
    }

    fn approx_eq_number(&self, other: &Self, tolerance: &Self) -> bool {
        for (i, c) in self.coeffs.iter().enumerate().take(INLINE_COEFF_COUNT) {
            if !c.approx_eq_number(&other.coeffs[i], &tolerance.coeffs[0]) {
                return false;
            }
        }
        true
    }

    fn total_cmp(&self, other: &Self) -> Ordering {
        for (i, c) in self.coeffs.iter().enumerate().take(INLINE_COEFF_COUNT) {
            let cmp = c.total_cmp(&other.coeffs[i]);
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        Ordering::Equal
    }
}
