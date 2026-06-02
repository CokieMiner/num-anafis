use crate::error::NumAnafisError;
use crate::number::{r as raw_r, s as raw_s};
use crate::{Number, Scalar};
use core::cmp::Ordering;
use core::hash::{Hash, Hasher};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::hash_map::DefaultHasher;

#[cfg(feature = "clifford")]
use alloc::vec::Vec;

#[cfg(feature = "clifford")]
use crate::clifford::{
    CliffordNumber, GeneratorSet, cga_gens, ci, e_minus, e_plus, e1, e2, e3, eps, inf, orig,
    pseudo3d, pseudo5d, qi, qj, qk, sj,
};

#[cfg(feature = "clifford")]
fn map_error(e: &NumAnafisError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

// ============================================================================
// PyScalar
// ============================================================================

#[pyclass(name = "Scalar", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyScalar(Scalar);

#[pymethods]
impl PyScalar {
    #[new]
    fn new(x: f64) -> Self {
        Self(raw_s(x))
    }

    #[staticmethod]
    fn from_ints(num: i64, den: i64) -> Self {
        Self(raw_r(num, den))
    }

    #[staticmethod]
    fn set_precision(bits: u32) -> bool {
        Scalar::set_precision(bits)
    }

    #[staticmethod]
    fn get_precision() -> u32 {
        Scalar::get_precision()
    }

    #[staticmethod]
    fn epsilon() -> Self {
        Self(Scalar::epsilon())
    }

    fn __repr__(&self) -> String {
        format!("Scalar({})", self.0)
    }

    fn __str__(&self) -> String {
        format!("{}", self.0)
    }

    fn __int__(&self, py: Python) -> PyResult<Py<PyAny>> {
        let s = match self.0.to_int() {
            Some(v) => format!("{v}"),
            None => return Err(PyValueError::new_err("Scalar is not an exact integer")),
        };
        let builtins = py.import("builtins")?;
        builtins.getattr("int")?.call1((s,)).map(Bound::unbind)
    }

    fn __float__(&self) -> f64 {
        let s = format!("{}", self.0);
        s.parse().unwrap_or(f64::NAN)
    }

    fn __bool__(&self) -> bool {
        !self.0.is_zero()
    }

    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.0.hash(&mut hasher);
        hasher.finish()
    }

    fn to_int(&self, py: Python) -> Option<Py<PyAny>> {
        let s = format!("{}", self.0.to_int()?);
        let builtins = py.import("builtins").ok()?;
        builtins
            .getattr("int")
            .ok()?
            .call1((s,))
            .ok()
            .map(Bound::unbind)
    }

    // --- arithmetic ---

    fn __add__(&self, other: &Self) -> Self {
        Self(&self.0 + &other.0)
    }
    fn __sub__(&self, other: &Self) -> Self {
        Self(&self.0 - &other.0)
    }
    fn __mul__(&self, other: &Self) -> Self {
        Self(&self.0 * &other.0)
    }
    fn __truediv__(&self, other: &Self) -> Self {
        Self(&self.0 / &other.0)
    }
    fn __neg__(&self) -> Self {
        Self(-&self.0)
    }
    fn __abs__(&self) -> Self {
        Self(self.0.abs())
    }
    fn __pow__(&self, exp: &Self, _mod: Option<&Self>) -> Self {
        Self(self.0.pow(&exp.0))
    }
    fn gcd(&self, other: &Self) -> Self {
        Self(self.0.gcd(&other.0))
    }

    // --- comparisons ---

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }
    fn __ne__(&self, other: &Self) -> bool {
        self.0 != other.0
    }
    fn __lt__(&self, other: &Self) -> bool {
        self.0.partial_cmp(&other.0) == Some(Ordering::Less)
    }
    fn __le__(&self, other: &Self) -> bool {
        matches!(
            self.0.partial_cmp(&other.0),
            Some(Ordering::Less | Ordering::Equal)
        )
    }
    fn __gt__(&self, other: &Self) -> bool {
        self.0.partial_cmp(&other.0) == Some(Ordering::Greater)
    }
    fn __ge__(&self, other: &Self) -> bool {
        matches!(
            self.0.partial_cmp(&other.0),
            Some(Ordering::Greater | Ordering::Equal)
        )
    }

    // --- Scalar-specific ---

    fn is_even(&self) -> bool {
        self.0.is_even()
    }
    fn is_odd(&self) -> bool {
        self.0.is_odd()
    }
    fn max(&self, other: &Self) -> Self {
        Self(self.0.clone().max(other.0.clone()))
    }
    fn min(&self, other: &Self) -> Self {
        Self(self.0.clone().min(other.0.clone()))
    }
    fn clamp(&self, lo: &Self, hi: &Self) -> Self {
        Self(self.0.clone().clamp(lo.0.clone(), hi.0.clone()))
    }

    // === Number trait methods ===

    // --- Trigonometric ---
    fn sin(&self) -> Self {
        Self(self.0.sin())
    }
    fn cos(&self) -> Self {
        Self(self.0.cos())
    }
    fn tan(&self) -> Self {
        Self(self.0.tan())
    }
    fn cot(&self) -> Self {
        Self(self.0.cot())
    }
    fn sec(&self) -> Self {
        Self(self.0.sec())
    }
    fn csc(&self) -> Self {
        Self(self.0.csc())
    }

    // --- Inverse Trigonometric ---
    fn asin(&self) -> Self {
        Self(self.0.asin())
    }
    fn acos(&self) -> Self {
        Self(self.0.acos())
    }
    fn atan(&self) -> Self {
        Self(self.0.atan())
    }
    fn acot(&self) -> Self {
        Self(self.0.acot())
    }
    fn asec(&self) -> Self {
        Self(self.0.asec())
    }
    fn acsc(&self) -> Self {
        Self(self.0.acsc())
    }

    // --- Hyperbolic ---
    fn sinh(&self) -> Self {
        Self(self.0.sinh())
    }
    fn cosh(&self) -> Self {
        Self(self.0.cosh())
    }
    fn tanh(&self) -> Self {
        Self(self.0.tanh())
    }
    fn coth(&self) -> Self {
        Self(self.0.coth())
    }
    fn sech(&self) -> Self {
        Self(self.0.sech())
    }
    fn csch(&self) -> Self {
        Self(self.0.csch())
    }

    // --- Inverse Hyperbolic ---
    fn asinh(&self) -> Self {
        Self(self.0.asinh())
    }
    fn acosh(&self) -> Self {
        Self(self.0.acosh())
    }
    fn atanh(&self) -> Self {
        Self(self.0.atanh())
    }
    fn acoth(&self) -> Self {
        Self(self.0.acoth())
    }
    fn asech(&self) -> Self {
        Self(self.0.asech())
    }
    fn acsch(&self) -> Self {
        Self(self.0.acsch())
    }

    // --- Exponential & Logarithmic ---
    fn exp(&self) -> Self {
        Self(self.0.exp())
    }
    fn expm1(&self) -> Self {
        Self(self.0.expm1())
    }
    fn exp_neg(&self) -> Self {
        Self(self.0.exp_neg())
    }
    fn ln(&self) -> Self {
        Self(self.0.ln())
    }
    fn log1p(&self) -> Self {
        Self(self.0.log1p())
    }

    // --- Powers / Roots ---
    fn sqrt(&self) -> Self {
        Self(self.0.sqrt())
    }
    fn cbrt(&self) -> Self {
        Self(self.0.cbrt())
    }

    // --- Basic Math ---
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    fn signum(&self) -> Self {
        Self(self.0.signum())
    }
    fn floor(&self) -> Self {
        Self(self.0.floor())
    }
    fn ceil(&self) -> Self {
        Self(self.0.ceil())
    }
    fn round(&self) -> Self {
        Self(self.0.round())
    }
    fn fract(&self) -> Self {
        Self(self.0.fract())
    }
    fn negate(&self) -> Self {
        Self(self.0.negate())
    }

    // --- Special Functions ---
    fn erf(&self) -> Self {
        Self(self.0.erf())
    }
    fn erfc(&self) -> Self {
        Self(self.0.erfc())
    }
    fn gamma(&self) -> Self {
        Self(self.0.gamma())
    }
    fn lgamma(&self) -> Self {
        Self(self.0.lgamma())
    }
    fn digamma(&self) -> Self {
        Self(self.0.digamma())
    }
    fn trigamma(&self) -> Self {
        Self(self.0.trigamma())
    }
    fn tetragamma(&self) -> Self {
        Self(self.0.tetragamma())
    }
    fn sinc(&self) -> Self {
        Self(self.0.sinc())
    }
    fn elliptic_k(&self) -> Self {
        Self(self.0.elliptic_k())
    }
    fn elliptic_e(&self) -> Self {
        Self(self.0.elliptic_e())
    }
    fn zeta(&self) -> Self {
        Self(self.0.zeta())
    }
    fn exp_polar(&self) -> Self {
        Self(self.0.exp_polar())
    }

    // --- Binary ---
    fn atan2(&self, x: &Self) -> Self {
        Self(self.0.atan2(&x.0))
    }
    fn log_base(&self, base: &Self) -> Self {
        Self(self.0.log_base(&base.0))
    }
    fn pow(&self, exp: &Self) -> Self {
        Self(self.0.pow(&exp.0))
    }
    fn besselj(&self, n: &Self) -> Self {
        Self(self.0.besselj(&n.0))
    }
    fn bessely(&self, n: &Self) -> Self {
        Self(self.0.bessely(&n.0))
    }
    fn besseli(&self, n: &Self) -> Self {
        Self(self.0.besseli(&n.0))
    }
    fn besselk(&self, n: &Self) -> Self {
        Self(self.0.besselk(&n.0))
    }
    fn polygamma(&self, n: &Self) -> Self {
        Self(self.0.polygamma(&n.0))
    }
    fn beta(&self, other: &Self) -> Self {
        Self(self.0.beta(&other.0))
    }
    fn zeta_deriv(&self, n: &Self) -> Self {
        Self(self.0.zeta_deriv(&n.0))
    }
    fn lambertw(&self, n: &Self) -> Self {
        Self(self.0.lambertw(&n.0))
    }
    fn hermite(&self, n: &Self) -> Self {
        Self(self.0.hermite(&n.0))
    }
    fn assoc_legendre(&self, l: &Self, m: &Self) -> Self {
        Self(self.0.assoc_legendre(&l.0, &m.0))
    }
    fn spherical_harmonic(&self, l: &Self, m: &Self, phi: &Self) -> Self {
        Self(self.0.spherical_harmonic(&l.0, &m.0, &phi.0))
    }

    // --- Properties ---
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn is_one(&self) -> bool {
        self.0.is_one()
    }
    fn is_neg_one(&self) -> bool {
        self.0.is_neg_one()
    }
    fn is_integer(&self) -> bool {
        self.0.is_integer()
    }
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }
    fn is_negative(&self) -> bool {
        self.0.is_negative()
    }
    fn is_positive(&self) -> bool {
        self.0.is_positive()
    }
    fn to_float(&self) -> Self {
        Self(self.0.to_float())
    }
    fn approx_eq(&self, other: &Self, tol: &Self) -> bool {
        self.0.approx_eq_number(&other.0, &tol.0)
    }
    fn total_cmp(&self, other: &Self) -> i8 {
        match self.0.total_cmp(&other.0) {
            Ordering::Less => -1,
            Ordering::Equal => 0,
            Ordering::Greater => 1,
        }
    }
    fn num_max(&self, other: &Self) -> Self {
        Self(self.0.num_max(&other.0))
    }
    fn num_min(&self, other: &Self) -> Self {
        Self(self.0.num_min(&other.0))
    }
}

// ============================================================================
// PyCliffordNumber
// ============================================================================

#[cfg(feature = "clifford")]
#[pyclass(name = "CliffordNumber", skip_from_py_object)]
#[derive(Clone)]
struct PyCliffordNumber(CliffordNumber);

#[cfg(feature = "clifford")]
#[pymethods]
impl PyCliffordNumber {
    // --- Constructors ---

    #[staticmethod]
    fn zero(gens: &PyGeneratorSet) -> PyResult<Self> {
        CliffordNumber::zero(gens.0.clone())
            .map(Self)
            .map_err(|e| map_error(&e))
    }

    #[staticmethod]
    fn scalar(gens: &PyGeneratorSet, value: &PyScalar) -> PyResult<Self> {
        CliffordNumber::scalar(gens.0.clone(), value.0.clone())
            .map(Self)
            .map_err(|e| map_error(&e))
    }

    #[staticmethod]
    fn generator(gens: &PyGeneratorSet, index: u8) -> PyResult<Self> {
        CliffordNumber::generator(gens.0.clone(), index)
            .map(Self)
            .map_err(|e| map_error(&e))
    }

    // --- Display ---

    fn __repr__(&self) -> String {
        format!("CliffordNumber({})", self.0)
    }

    fn __str__(&self) -> String {
        format!("{}", self.0)
    }

    // --- Accessors ---

    fn blade_count(&self) -> usize {
        self.0.blade_count()
    }
    const fn n_generators(&self) -> usize {
        self.0.n_generators()
    }
    fn generator_set(&self) -> PyGeneratorSet {
        PyGeneratorSet(self.0.generator_set().clone())
    }
    fn coeff(&self, blade: usize) -> PyScalar {
        PyScalar(self.0.coeff(blade).clone())
    }
    fn set_coeff(&mut self, blade: usize, value: &PyScalar) {
        self.0.set_coeff(blade, value.0.clone());
    }
    fn nonzero_blades(&self) -> Vec<(usize, PyScalar)> {
        self.0
            .nonzero_blades()
            .map(|(b, c)| (b, PyScalar(c.clone())))
            .collect()
    }
    const fn is_heap_allocated(&self) -> bool {
        self.0.is_heap_allocated()
    }
    fn chop(&self) -> Self {
        Self(self.0.clone().chop())
    }

    // --- Geometric operations ---

    fn geometric_mul(&self, other: &Self) -> Self {
        Self(self.0.geometric_mul(&other.0))
    }
    fn outer_product(&self, other: &Self) -> Self {
        Self(self.0.outer_product(&other.0))
    }
    fn inner_product(&self, other: &Self) -> Self {
        Self(self.0.inner_product(&other.0))
    }
    fn scalar_product(&self, other: &Self) -> PyScalar {
        PyScalar(self.0.scalar_product(&other.0))
    }
    fn geometric_inverse(&self) -> Self {
        Self(self.0.geometric_inverse())
    }

    // --- Grade operations ---

    fn grade(&self, k: u32) -> Self {
        Self(self.0.grade(k))
    }
    fn reverse(&self) -> Self {
        Self(self.0.reverse())
    }
    fn grade_involution(&self) -> Self {
        Self(self.0.grade_involution())
    }
    fn clifford_conjugate(&self) -> Self {
        Self(self.0.clifford_conjugate())
    }
    fn norm_sq(&self) -> PyScalar {
        PyScalar(self.0.norm_sq())
    }

    // --- Arithmetic ---

    fn __add__(&self, other: &Self) -> Self {
        Self(&self.0 + &other.0)
    }
    fn __sub__(&self, other: &Self) -> Self {
        Self(&self.0 - &other.0)
    }
    fn __mul__(&self, other: &Self) -> Self {
        Self(self.0.geometric_mul(&other.0))
    }
    fn __truediv__(&self, other: &Self) -> Self {
        Self(&self.0 / &other.0)
    }
    fn __neg__(&self) -> Self {
        Self(-&self.0)
    }
    fn __abs__(&self) -> Self {
        Self(self.0.abs())
    }
    fn __pow__(&self, exp: &Self, _mod: Option<&Self>) -> Self {
        Self(self.0.pow(&exp.0))
    }

    // --- Comparisons ---

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }
    fn __ne__(&self, other: &Self) -> bool {
        self.0 != other.0
    }

    // === Number trait methods ===

    // --- Trigonometric ---
    fn sin(&self) -> Self {
        Self(self.0.sin())
    }
    fn cos(&self) -> Self {
        Self(self.0.cos())
    }
    fn tan(&self) -> Self {
        Self(self.0.tan())
    }
    fn cot(&self) -> Self {
        Self(self.0.cot())
    }
    fn sec(&self) -> Self {
        Self(self.0.sec())
    }
    fn csc(&self) -> Self {
        Self(self.0.csc())
    }

    // --- Inverse Trigonometric ---
    fn asin(&self) -> Self {
        Self(self.0.asin())
    }
    fn acos(&self) -> Self {
        Self(self.0.acos())
    }
    fn atan(&self) -> Self {
        Self(self.0.atan())
    }
    fn acot(&self) -> Self {
        Self(self.0.acot())
    }
    fn asec(&self) -> Self {
        Self(self.0.asec())
    }
    fn acsc(&self) -> Self {
        Self(self.0.acsc())
    }

    // --- Hyperbolic ---
    fn sinh(&self) -> Self {
        Self(self.0.sinh())
    }
    fn cosh(&self) -> Self {
        Self(self.0.cosh())
    }
    fn tanh(&self) -> Self {
        Self(self.0.tanh())
    }
    fn coth(&self) -> Self {
        Self(self.0.coth())
    }
    fn sech(&self) -> Self {
        Self(self.0.sech())
    }
    fn csch(&self) -> Self {
        Self(self.0.csch())
    }

    // --- Inverse Hyperbolic ---
    fn asinh(&self) -> Self {
        Self(self.0.asinh())
    }
    fn acosh(&self) -> Self {
        Self(self.0.acosh())
    }
    fn atanh(&self) -> Self {
        Self(self.0.atanh())
    }
    fn acoth(&self) -> Self {
        Self(self.0.acoth())
    }
    fn asech(&self) -> Self {
        Self(self.0.asech())
    }
    fn acsch(&self) -> Self {
        Self(self.0.acsch())
    }

    // --- Exponential & Logarithmic ---
    fn exp(&self) -> Self {
        Self(self.0.exp())
    }
    fn expm1(&self) -> Self {
        Self(self.0.expm1())
    }
    fn exp_neg(&self) -> Self {
        Self(self.0.exp_neg())
    }
    fn ln(&self) -> Self {
        Self(self.0.ln())
    }
    fn log1p(&self) -> Self {
        Self(self.0.log1p())
    }

    // --- Powers / Roots ---
    fn sqrt(&self) -> Self {
        Self(self.0.sqrt())
    }
    fn cbrt(&self) -> Self {
        Self(self.0.cbrt())
    }

    // --- Basic Math ---
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    fn signum(&self) -> Self {
        Self(self.0.signum())
    }
    fn floor(&self) -> Self {
        Self(self.0.floor())
    }
    fn ceil(&self) -> Self {
        Self(self.0.ceil())
    }
    fn round(&self) -> Self {
        Self(self.0.round())
    }
    fn fract(&self) -> Self {
        Self(self.0.fract())
    }
    fn negate(&self) -> Self {
        Self(self.0.negate())
    }

    // --- Special Functions ---
    fn erf(&self) -> Self {
        Self(self.0.erf())
    }
    fn erfc(&self) -> Self {
        Self(self.0.erfc())
    }
    fn gamma(&self) -> Self {
        Self(self.0.gamma())
    }
    fn lgamma(&self) -> Self {
        Self(self.0.lgamma())
    }
    fn digamma(&self) -> Self {
        Self(self.0.digamma())
    }
    fn trigamma(&self) -> Self {
        Self(self.0.trigamma())
    }
    fn tetragamma(&self) -> Self {
        Self(self.0.tetragamma())
    }
    fn sinc(&self) -> Self {
        Self(self.0.sinc())
    }
    fn elliptic_k(&self) -> Self {
        Self(self.0.elliptic_k())
    }
    fn elliptic_e(&self) -> Self {
        Self(self.0.elliptic_e())
    }
    fn zeta(&self) -> Self {
        Self(self.0.zeta())
    }
    fn exp_polar(&self) -> Self {
        Self(self.0.exp_polar())
    }

    // --- Binary ---
    fn atan2(&self, x: &Self) -> Self {
        Self(self.0.atan2(&x.0))
    }
    fn log_base(&self, base: &Self) -> Self {
        Self(self.0.log_base(&base.0))
    }
    fn pow(&self, exp: &Self) -> Self {
        Self(self.0.pow(&exp.0))
    }
    fn besselj(&self, n: &Self) -> Self {
        Self(self.0.besselj(&n.0))
    }
    fn bessely(&self, n: &Self) -> Self {
        Self(self.0.bessely(&n.0))
    }
    fn besseli(&self, n: &Self) -> Self {
        Self(self.0.besseli(&n.0))
    }
    fn besselk(&self, n: &Self) -> Self {
        Self(self.0.besselk(&n.0))
    }
    fn polygamma(&self, n: &Self) -> Self {
        Self(self.0.polygamma(&n.0))
    }
    fn beta(&self, other: &Self) -> Self {
        Self(self.0.beta(&other.0))
    }
    fn zeta_deriv(&self, n: &Self) -> Self {
        Self(self.0.zeta_deriv(&n.0))
    }
    fn lambertw(&self, n: &Self) -> Self {
        Self(self.0.lambertw(&n.0))
    }
    fn hermite(&self, n: &Self) -> Self {
        Self(self.0.hermite(&n.0))
    }
    fn assoc_legendre(&self, l: &Self, m: &Self) -> Self {
        Self(self.0.assoc_legendre(&l.0, &m.0))
    }
    fn spherical_harmonic(&self, l: &Self, m: &Self, phi: &Self) -> Self {
        Self(self.0.spherical_harmonic(&l.0, &m.0, &phi.0))
    }

    // --- Properties ---
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn is_one(&self) -> bool {
        self.0.is_one()
    }
    fn is_neg_one(&self) -> bool {
        self.0.is_neg_one()
    }
    fn is_integer(&self) -> bool {
        self.0.is_integer()
    }
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }
    fn is_negative(&self) -> bool {
        self.0.is_negative()
    }
    fn is_positive(&self) -> bool {
        self.0.is_positive()
    }
    fn to_float(&self) -> Self {
        Self(self.0.to_float())
    }
    fn approx_eq(&self, other: &Self, tol: &Self) -> bool {
        self.0.approx_eq_number(&other.0, &tol.0)
    }
    fn total_cmp(&self, other: &Self) -> i8 {
        match self.0.total_cmp(&other.0) {
            Ordering::Less => -1,
            Ordering::Equal => 0,
            Ordering::Greater => 1,
        }
    }
    fn num_max(&self, other: &Self) -> Self {
        Self(self.0.num_max(&other.0))
    }
    fn num_min(&self, other: &Self) -> Self {
        Self(self.0.num_min(&other.0))
    }
}

// ============================================================================
// PyGeneratorSet
// ============================================================================

#[cfg(feature = "clifford")]
#[pyclass(name = "GeneratorSet", skip_from_py_object)]
#[derive(Clone)]
struct PyGeneratorSet(GeneratorSet);

#[cfg(feature = "clifford")]
#[pymethods]
impl PyGeneratorSet {
    #[staticmethod]
    fn empty() -> Self {
        Self(GeneratorSet::empty())
    }

    #[staticmethod]
    fn cga() -> Self {
        Self(cga_gens())
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.0)
    }
    fn __str__(&self) -> String {
        format!("{:?}", self.0)
    }
    const fn __len__(&self) -> usize {
        self.0.len()
    }
    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }
    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.0.hash(&mut hasher);
        hasher.finish()
    }

    const fn len(&self) -> usize {
        self.0.len()
    }
    const fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    fn id_at(&self, i: usize) -> u32 {
        self.0.id_at(i)
    }
    fn metric_at(&self, i: usize) -> i8 {
        self.0.metric_at(i)
    }
}

// ============================================================================
// Clifford free functions
// ============================================================================

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e1() -> PyCliffordNumber {
    PyCliffordNumber(e1())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e2() -> PyCliffordNumber {
    PyCliffordNumber(e2())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e3() -> PyCliffordNumber {
    PyCliffordNumber(e3())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e_plus() -> PyCliffordNumber {
    PyCliffordNumber(e_plus())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e_minus() -> PyCliffordNumber {
    PyCliffordNumber(e_minus())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_orig() -> PyCliffordNumber {
    PyCliffordNumber(orig())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_inf() -> PyCliffordNumber {
    PyCliffordNumber(inf())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_qi() -> PyCliffordNumber {
    PyCliffordNumber(qi())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_qj() -> PyCliffordNumber {
    PyCliffordNumber(qj())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_qk() -> PyCliffordNumber {
    PyCliffordNumber(qk())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_pseudo3d() -> PyCliffordNumber {
    PyCliffordNumber(pseudo3d())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_pseudo5d() -> PyCliffordNumber {
    PyCliffordNumber(pseudo5d())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_eps() -> PyCliffordNumber {
    PyCliffordNumber(eps())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_ci() -> PyCliffordNumber {
    PyCliffordNumber(ci())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_sj() -> PyCliffordNumber {
    PyCliffordNumber(sj())
}

// ============================================================================
// Common free functions
// ============================================================================

#[pyfunction]
fn s(x: f64) -> PyScalar {
    PyScalar(raw_s(x))
}

#[pyfunction]
fn r(num: i64, den: i64) -> PyScalar {
    PyScalar(raw_r(num, den))
}

#[pyfunction]
fn get_precision() -> u32 {
    Scalar::get_precision()
}

// ============================================================================
// Module registration
// ============================================================================

#[pymodule]
fn num_anafis_py(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyScalar>()?;
    m.add_function(wrap_pyfunction!(s, m)?)?;
    m.add_function(wrap_pyfunction!(r, m)?)?;
    m.add_function(wrap_pyfunction!(get_precision, m)?)?;

    #[cfg(feature = "clifford")]
    {
        m.add_class::<PyCliffordNumber>()?;
        m.add_class::<PyGeneratorSet>()?;
        m.add_function(wrap_pyfunction!(py_e1, m)?)?;
        m.add_function(wrap_pyfunction!(py_e2, m)?)?;
        m.add_function(wrap_pyfunction!(py_e3, m)?)?;
        m.add_function(wrap_pyfunction!(py_e_plus, m)?)?;
        m.add_function(wrap_pyfunction!(py_e_minus, m)?)?;
        m.add_function(wrap_pyfunction!(py_orig, m)?)?;
        m.add_function(wrap_pyfunction!(py_inf, m)?)?;
        m.add_function(wrap_pyfunction!(py_qi, m)?)?;
        m.add_function(wrap_pyfunction!(py_qj, m)?)?;
        m.add_function(wrap_pyfunction!(py_qk, m)?)?;
        m.add_function(wrap_pyfunction!(py_pseudo3d, m)?)?;
        m.add_function(wrap_pyfunction!(py_pseudo5d, m)?)?;
        m.add_function(wrap_pyfunction!(py_eps, m)?)?;
        m.add_function(wrap_pyfunction!(py_ci, m)?)?;
        m.add_function(wrap_pyfunction!(py_sj, m)?)?;
    }

    Ok(())
}
