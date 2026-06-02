use alloc::boxed::Box;
use alloc::vec;
use alloc::vec::Vec;

use super::CliffordNumber;
use super::large_mat::{
    apply_via_spectral as large_apply_via_spectral,
    apply_via_spectral_closure as large_apply_via_spectral_closure,
};
use super::matrix::{
    Cmplx, Mat2C, can_embed, compute_blade_basis, eigendecompose, from_matrix_with_gens,
    to_matrix_with_basis,
};
use super::types::{GeneratorSet, scalar_nan, scalar_pi, scalar_zero};
use crate::number::s;
use crate::scalar::Scalar;
use crate::traits::Number;
use core::cmp::Ordering;
use core::ops::Add;

#[derive(Debug, Clone, PartialEq)]
pub enum SpectralFn {
    // Trigonometric
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Asin,
    Acos,
    Atan,
    Acot,
    Asec,
    Acsc,

    // Hyperbolic
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Sech,
    Csch,
    Asinh,
    Acosh,
    Atanh,
    Acoth,
    Asech,
    Acsch,

    // Exponential & Power
    Exp,
    Expm1,
    ExpNeg,
    Ln,
    Log1p,
    Sqrt,
    Cbrt,

    // Special
    Erf,
    Erfc,
    Gamma,
    LGamma,
    Polygamma(usize),
    Zeta,
    ZetaDeriv(usize),
    EllipticK,
    EllipticE,
    Sinc,

    // AST Nodes
    Var,
    Const(Cmplx),
    Add(Vec<Self>),
    Mul(Vec<Self>),
    Inv(Box<Self>),
    Neg(Box<Self>),
    Compose(Box<Self>, Box<Self>),
}

impl SpectralFn {
    #[allow(
        clippy::pattern_type_mismatch,
        reason = "match ergonomics produce cleaner code here"
    )]
    pub fn eval(&self, c: &Cmplx) -> Cmplx {
        match self {
            Self::Sin => c.sin(),
            Self::Cos => c.cos(),
            Self::Tan => c.sin().div(&c.cos()),
            Self::Cot => c.cos().div(&c.sin()),
            Self::Sec => Cmplx::one().div(&c.cos()),
            Self::Csc => Cmplx::one().div(&c.sin()),
            Self::Asin => c.apply_real_only(&|s| s.asin()),
            Self::Acos => c.apply_real_only(&|s| s.acos()),
            Self::Atan => c.apply_real_only(&|s| s.atan()),
            Self::Acot => c.apply_real_only(&|s| s.acot()),
            Self::Asec => c.apply_real_only(&|s| s.asec()),
            Self::Acsc => c.apply_real_only(&|s| s.acsc()),

            Self::Sinh => c.sinh(),
            Self::Cosh => c.cosh(),
            Self::Tanh => c.sinh().div(&c.cosh()),
            Self::Coth => c.cosh().div(&c.sinh()),
            Self::Sech => Cmplx::one().div(&c.cosh()),
            Self::Csch => Cmplx::one().div(&c.sinh()),
            Self::Asinh => c.apply_real_only(&|s| s.asinh()),
            Self::Acosh => c.apply_real_only(&|s| s.acosh()),
            Self::Atanh => c.apply_real_only(&|s| s.atanh()),
            Self::Acoth => c.apply_real_only(&|s| s.acoth()),
            Self::Asech => c.apply_real_only(&|s| s.asech()),
            Self::Acsch => c.apply_real_only(&|s| s.acsch()),

            Self::Exp => c.exp(),
            Self::Expm1 => c.exp().sub(&Cmplx::one()),
            Self::ExpNeg => c.neg().exp(),
            Self::Ln => c.ln(),
            Self::Log1p => c.add(&Cmplx::one()).ln(),
            Self::Sqrt => c.sqrt(),
            Self::Cbrt => c.apply_real_only(&|s| s.cbrt()),

            Self::Erf => c.apply_real_only(&|s| s.erf()),
            Self::Erfc => c.apply_real_only(&|s| s.erfc()),
            Self::Gamma => c.apply_real_only(&|s| s.gamma()),
            Self::LGamma => c.apply_real_only(&|s| s.lgamma()),
            Self::Polygamma(n) => {
                let n_scalar = s(*n);
                c.apply_real_only(&|s| s.polygamma(&n_scalar))
            }
            Self::Zeta => c.apply_real_only(&|s| s.zeta()),
            Self::ZetaDeriv(n) => {
                let n_scalar = s(*n);
                c.apply_real_only(&|s| s.zeta_deriv(&n_scalar))
            }
            Self::EllipticK => c.apply_real_only(&|s| s.elliptic_k()),
            Self::EllipticE => c.apply_real_only(&|s| s.elliptic_e()),
            Self::Sinc => {
                if c.0.is_zero() && c.1.is_zero() {
                    Cmplx::one()
                } else {
                    c.sin().div(c)
                }
            }

            Self::Var => c.clone(),
            Self::Const(v) => v.clone(),
            Self::Add(args) => {
                let mut sum = Cmplx::zero();
                for a in args {
                    sum = sum.add(&a.eval(c));
                }
                sum
            }
            Self::Mul(args) => {
                let mut prod = Cmplx::one();
                for a in args {
                    prod = prod.mul(&a.eval(c));
                }
                prod
            }
            Self::Inv(a) => Cmplx::one().div(&a.eval(c)),
            Self::Neg(a) => a.eval(c).neg(),
            Self::Compose(outer, inner) => outer.eval(&inner.eval(c)),
        }
    }

    #[allow(
        clippy::too_many_lines,
        clippy::pattern_type_mismatch,
        reason = "Large match block and match ergonomics are much cleaner here"
    )]
    pub fn derivative(&self) -> Self {
        match self {
            Self::Sin => Self::Cos,
            Self::Cos => Self::Neg(Box::new(Self::Sin)),
            Self::Tan => Self::Mul(vec![Self::Sec, Self::Sec]),
            Self::Cot => Self::Neg(Box::new(Self::Mul(vec![Self::Csc, Self::Csc]))),
            Self::Sec => Self::Mul(vec![Self::Sec, Self::Tan]),
            Self::Csc => Self::Neg(Box::new(Self::Mul(vec![Self::Csc, Self::Cot]))),

            Self::Asin => Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Const(Cmplx::one()),
                    Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
                ])),
            ))),
            Self::Acos => Self::Neg(Box::new(Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Const(Cmplx::one()),
                    Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
                ])),
            ))))),
            Self::Atan => Self::Inv(Box::new(Self::Add(vec![
                Self::Const(Cmplx::one()),
                Self::Mul(vec![Self::Var, Self::Var]),
            ]))),
            Self::Acot => Self::Neg(Box::new(Self::Inv(Box::new(Self::Add(vec![
                Self::Const(Cmplx::one()),
                Self::Mul(vec![Self::Var, Self::Var]),
            ]))))),
            Self::Asec => Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Mul(vec![Self::Var, Self::Var]),
                        Self::Const(Cmplx::one().neg()),
                    ])),
                ),
            ]))),
            Self::Acsc => Self::Neg(Box::new(Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Mul(vec![Self::Var, Self::Var]),
                        Self::Const(Cmplx::one().neg()),
                    ])),
                ),
            ]))))),

            Self::Sinh => Self::Cosh,
            Self::Cosh => Self::Sinh,
            Self::Tanh => Self::Mul(vec![Self::Sech, Self::Sech]),
            Self::Coth => Self::Neg(Box::new(Self::Mul(vec![Self::Csch, Self::Csch]))),
            Self::Sech => Self::Neg(Box::new(Self::Mul(vec![Self::Sech, Self::Tanh]))),
            Self::Csch => Self::Neg(Box::new(Self::Mul(vec![Self::Csch, Self::Coth]))),

            Self::Asinh => Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Mul(vec![Self::Var, Self::Var]),
                    Self::Const(Cmplx::one()),
                ])),
            ))),
            Self::Acosh => Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Mul(vec![Self::Var, Self::Var]),
                    Self::Const(Cmplx::one().neg()),
                ])),
            ))),
            Self::Atanh | Self::Acoth => Self::Inv(Box::new(Self::Add(vec![
                Self::Const(Cmplx::one()),
                Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
            ]))),
            Self::Asech => Self::Neg(Box::new(Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Const(Cmplx::one()),
                        Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
                    ])),
                ),
            ]))))),
            Self::Acsch => Self::Neg(Box::new(Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Const(Cmplx::one()),
                        Self::Mul(vec![Self::Var, Self::Var]),
                    ])),
                ),
            ]))))),

            Self::Exp | Self::Expm1 => Self::Exp,
            Self::ExpNeg => Self::Neg(Box::new(Self::ExpNeg)),
            Self::Ln => Self::Inv(Box::new(Self::Var)),
            Self::Log1p => Self::Inv(Box::new(Self::Add(vec![
                Self::Var,
                Self::Const(Cmplx::one()),
            ]))),

            Self::Sqrt => Self::Mul(vec![
                Self::Const(Cmplx::new(s(0.5), scalar_zero())),
                Self::Inv(Box::new(Self::Sqrt)),
            ]),
            Self::Cbrt => Self::Mul(vec![
                Self::Const(Cmplx::new(s(1.0 / 3.0), scalar_zero())),
                Self::Inv(Box::new(Self::Mul(vec![Self::Cbrt, Self::Cbrt]))),
            ]),

            Self::LGamma => Self::Polygamma(0),
            Self::Polygamma(n) => Self::Polygamma(*n + 1),
            Self::Gamma => Self::Mul(vec![Self::Gamma, Self::Polygamma(0)]),
            Self::Zeta => Self::ZetaDeriv(1),
            Self::ZetaDeriv(n) => Self::ZetaDeriv(*n + 1),

            Self::Sinc => Self::Mul(vec![
                Self::Inv(Box::new(Self::Var)),
                Self::Add(vec![Self::Cos, Self::Neg(Box::new(Self::Sinc))]),
            ]),

            Self::Erf => {
                // erf'(x) = (2/√π)·exp(-x²)
                let two_over_sqrt_pi = Cmplx::new(s(2) / scalar_pi().sqrt(), scalar_zero());
                Self::Mul(vec![
                    Self::Const(two_over_sqrt_pi),
                    Self::Compose(
                        Box::new(Self::ExpNeg),
                        Box::new(Self::Mul(vec![Self::Var, Self::Var])),
                    ),
                ])
            }
            Self::Erfc => {
                // erfc'(x) = -erf'(x)
                let two_over_sqrt_pi = Cmplx::new(s(2) / scalar_pi().sqrt(), scalar_zero());
                Self::Neg(Box::new(Self::Mul(vec![
                    Self::Const(two_over_sqrt_pi),
                    Self::Compose(
                        Box::new(Self::ExpNeg),
                        Box::new(Self::Mul(vec![Self::Var, Self::Var])),
                    ),
                ])))
            }
            Self::EllipticK => {
                // K'(k) = (E(k) − (1−k²)·K(k)) / (k·(1−k²))
                let one_m_k2 = Self::Add(vec![
                    Self::Const(Cmplx::one()),
                    Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
                ]);
                Self::Mul(vec![
                    Self::Inv(Box::new(Self::Var)),
                    Self::Inv(Box::new(one_m_k2.clone())),
                    Self::Add(vec![
                        Self::EllipticE,
                        Self::Neg(Box::new(Self::Mul(vec![one_m_k2, Self::EllipticK]))),
                    ]),
                ])
            }
            Self::EllipticE => {
                // E'(k) = (E(k) − K(k)) / k
                Self::Mul(vec![
                    Self::Inv(Box::new(Self::Var)),
                    Self::Add(vec![Self::EllipticE, Self::Neg(Box::new(Self::EllipticK))]),
                ])
            }

            // Base cases
            Self::Var => Self::Const(Cmplx::one()),
            Self::Const(_) => Self::Const(Cmplx::zero()),

            // AST Rules
            Self::Compose(outer, inner) => Self::Mul(vec![
                Self::Compose(Box::new(outer.derivative()), inner.clone()),
                inner.derivative(),
            ]),
            Self::Neg(a) => Self::Neg(Box::new(a.derivative())),
            Self::Inv(a) => Self::Neg(Box::new(Self::Mul(vec![
                a.derivative(),
                Self::Inv(Box::new(Self::Mul(vec![*a.clone(), *a.clone()]))),
            ]))),
            Self::Add(args) => Self::Add(args.iter().map(Self::derivative).collect()),
            Self::Mul(args) => {
                let mut sum_terms = Vec::with_capacity(args.len());
                for i in 0..args.len() {
                    let mut prod_terms = args.clone();
                    prod_terms[i] = prod_terms[i].derivative();
                    sum_terms.push(Self::Mul(prod_terms));
                }
                Self::Add(sum_terms)
            }
        }
    }
}

// ============================================================================
// Jet — truncated Taylor series (dual numbers of degree K)
//
// Each Jet represents a truncated Taylor series
//   f(x₀ + ε) = f(x₀) + f'(x₀) ε + f''(x₀)/2 ε² + … + f^(K)(x₀)/K! εᴷ
// where coeffs[i] = f^(i)(x₀) / i! .
//
// Evaluated once instead of calling `.derivative()` in a loop (which can
// grow the AST exponentially via the product/chain rules).
// ============================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Jet {
    pub coeffs: Vec<Cmplx>,
}

impl Jet {
    pub const fn new(coeffs: Vec<Cmplx>) -> Self {
        Self { coeffs }
    }

    pub const fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    pub fn constant(c: Cmplx, degree: usize) -> Self {
        let mut coeffs = vec![Cmplx::zero(); degree + 1];
        coeffs[0] = c;
        Self { coeffs }
    }

    pub fn zero(degree: usize) -> Self {
        Self::constant(Cmplx::zero(), degree)
    }

    pub fn one(degree: usize) -> Self {
        Self::constant(Cmplx::one(), degree)
    }

    /// Variable jet: `λ + ε` (used as the differentiation point in Block Parlett).
    pub fn variable(lambda: &Cmplx, degree: usize) -> Self {
        let mut coeffs = vec![Cmplx::zero(); degree + 1];
        coeffs[0] = lambda.clone();
        if degree >= 1 {
            coeffs[1] = Cmplx::one();
        }
        Self { coeffs }
    }

    // ------------------------------------------------------------------
    // Arithmetic
    // ------------------------------------------------------------------

    pub fn add(&self, other: &Self) -> Self {
        let max_deg = self.degree().max(other.degree());
        let mut coeffs = vec![Cmplx::zero(); max_deg + 1];
        for (i, c) in coeffs.iter_mut().enumerate().take(self.degree() + 1) {
            *c = c.add(&self.coeffs[i]);
        }
        for (i, c) in coeffs.iter_mut().enumerate().take(other.degree() + 1) {
            *c = c.add(&other.coeffs[i]);
        }
        Self { coeffs }
    }

    pub fn neg(&self) -> Self {
        Self {
            coeffs: self.coeffs.iter().map(Cmplx::neg).collect(),
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Scale every coefficient by a complex scalar.
    pub fn scale(&self, s: &Cmplx) -> Self {
        Self {
            coeffs: self.coeffs.iter().map(|c| c.mul(s)).collect(),
        }
    }

    /// Cauchy product of two truncated Taylor series:
    ///   (f·g)[n] = Σ_{i+j=n} f[i]·g[j]
    pub fn mul(&self, other: &Self) -> Self {
        let k = self.degree() + other.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        for i in 0..=self.degree() {
            for j in 0..=other.degree() {
                coeffs[i + j] = coeffs[i + j].add(&self.coeffs[i].mul(&other.coeffs[j]));
            }
        }
        Self { coeffs }
    }

    /// Multiplicative inverse via recurrence:
    ///   inv[0] = 1 / f[0]
    ///   inv[n] = -(1/f[0]) · Σ_{i=1}^{n} f[i] · inv[n-i]
    pub fn inv(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        let inv0 = Cmplx::one().div(&self.coeffs[0]);
        coeffs[0] = inv0.clone();
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for i in 1..=n {
                sum = sum.add(&self.coeffs[i].mul(&coeffs[n - i]));
            }
            coeffs[n] = sum.neg().mul(&inv0);
        }
        Self { coeffs }
    }

    /// Division: `a / b = a · b⁻¹`
    pub fn div(&self, other: &Self) -> Self {
        self.mul(&other.inv())
    }

    // ------------------------------------------------------------------
    // Transcendental recurrences
    // ------------------------------------------------------------------

    /// Exponential via recurrence:
    ///   exp[0] = exp(f[0])
    ///   exp[n] = (1/n) · Σ_{i=1}^{n} i · f[i] · exp[n-i]
    pub fn exp(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        coeffs[0] = self.coeffs[0].exp();
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for i in 1..=n {
                let i_c = Cmplx::new(s(i), scalar_zero());
                sum = sum.add(&i_c.mul(&self.coeffs[i]).mul(&coeffs[n - i]));
            }
            let n_c = Cmplx::new(s(n), scalar_zero());
            coeffs[n] = sum.div(&n_c);
        }
        Self { coeffs }
    }

    /// Natural log via recurrence:
    ///   ln[0] = ln(f[0])
    ///   ln[n] = (f[n] - (1/n) · Σ_{i=1}^{n-1} i · ln[i] · f[n-i]) / f[0]
    pub fn ln(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        coeffs[0] = self.coeffs[0].ln();
        let inv_f0 = Cmplx::one().div(&self.coeffs[0]);
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for (i, ci) in coeffs.iter().enumerate().take(n).skip(1) {
                let i_c = Cmplx::new(s(i), scalar_zero());
                sum = sum.add(&i_c.mul(ci).mul(&self.coeffs[n - i]));
            }
            coeffs[n] = self.coeffs[n]
                .sub(&sum.div(&Cmplx::new(s(n), scalar_zero())))
                .mul(&inv_f0);
        }
        Self { coeffs }
    }

    /// Sine (and simultaneously cosine) via mutual recurrence:
    ///   sin[0] = sin(f[0])   cos[0] = cos(f[0])
    ///   sin[n] = (1/n) · Σ i·f[i]·cos[n-i]
    ///   cos[n] = -(1/n) · Σ i·f[i]·sin[n-i]
    fn sin_cos(&self) -> (Self, Self) {
        let k = self.degree();
        let mut sin_c = vec![Cmplx::zero(); k + 1];
        let mut cos_c = vec![Cmplx::zero(); k + 1];
        sin_c[0] = self.coeffs[0].sin();
        cos_c[0] = self.coeffs[0].cos();
        for n in 1..=k {
            let mut s_sum = Cmplx::zero();
            let mut c_sum = Cmplx::zero();
            for i in 1..=n {
                let i_c = Cmplx::new(s(i), scalar_zero());
                let fi = self.coeffs[i].mul(&i_c);
                s_sum = s_sum.add(&fi.mul(&cos_c[n - i]));
                c_sum = c_sum.add(&fi.mul(&sin_c[n - i]));
            }
            let n_c = Cmplx::new(s(n), scalar_zero());
            sin_c[n] = s_sum.div(&n_c);
            cos_c[n] = c_sum.neg().div(&n_c);
        }
        (Self { coeffs: sin_c }, Self { coeffs: cos_c })
    }

    pub fn sin(&self) -> Self {
        self.sin_cos().0
    }

    pub fn cos(&self) -> Self {
        self.sin_cos().1
    }

    /// Hyperbolic sine (and simultaneously hyperbolic cosine) via mutual recurrence:
    ///   sinh[0] = sinh(f[0])  cosh[0] = cosh(f[0])
    ///   sinh[n] = (1/n) · Σ i·f[i]·cosh[n-i]
    ///   cosh[n] = (1/n) · Σ i·f[i]·sinh[n-i]
    fn sinh_cosh(&self) -> (Self, Self) {
        let k = self.degree();
        let mut sinh_c = vec![Cmplx::zero(); k + 1];
        let mut cosh_c = vec![Cmplx::zero(); k + 1];
        sinh_c[0] = self.coeffs[0].sinh();
        cosh_c[0] = self.coeffs[0].cosh();
        for n in 1..=k {
            let mut s_sum = Cmplx::zero();
            let mut c_sum = Cmplx::zero();
            for i in 1..=n {
                let i_c = Cmplx::new(s(i), scalar_zero());
                let fi = self.coeffs[i].mul(&i_c);
                s_sum = s_sum.add(&fi.mul(&cosh_c[n - i]));
                c_sum = c_sum.add(&fi.mul(&sinh_c[n - i]));
            }
            let n_c = Cmplx::new(s(n), scalar_zero());
            sinh_c[n] = s_sum.div(&n_c);
            cosh_c[n] = c_sum.div(&n_c);
        }
        (Self { coeffs: sinh_c }, Self { coeffs: cosh_c })
    }

    pub fn sinh(&self) -> Self {
        self.sinh_cosh().0
    }

    pub fn cosh(&self) -> Self {
        self.sinh_cosh().1
    }

    /// Square root via recurrence:
    ///   sqrt[0] = sqrt(f[0])
    ///   sqrt[n] = (f[n] - Σ_{i=1}^{n-1} sqrt[i]·sqrt[n-i]) / (2·sqrt[0])
    pub fn sqrt(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        coeffs[0] = self.coeffs[0].sqrt();
        let two_s0 = coeffs[0].add(&coeffs[0]); // 2·sqrt(f0)
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for i in 1..n {
                sum = sum.add(&coeffs[i].mul(&coeffs[n - i]));
            }
            coeffs[n] = self.coeffs[n].sub(&sum).div(&two_s0);
        }
        Self { coeffs }
    }
}

impl SpectralFn {
    /// Evaluate this AST over a **Jet** (truncated Taylor series).
    ///
    /// The result is a Jet whose coefficients are the Taylor coefficients of the
    /// function composed with the input series, obviating repeated calls to
    /// `.derivative()`.
    #[allow(
        clippy::pattern_type_mismatch,
        reason = "match ergonomics produce cleaner code here"
    )]
    pub fn eval_jet(&self, x: &Jet) -> Jet {
        let deg = x.degree();
        match self {
            // --- Primitives with direct Jet recurrences ---
            Self::Exp => x.exp(),
            Self::Expm1 => x.exp().sub(&Jet::one(deg)),
            Self::ExpNeg => x.neg().exp(),
            Self::Ln => x.ln(),
            Self::Log1p => {
                let one = Jet::one(deg);
                x.add(&one).ln()
            }
            Self::Sin => x.sin(),
            Self::Cos => x.cos(),
            Self::Tan => {
                let s = x.sin();
                let c = x.cos();
                s.div(&c)
            }
            Self::Cot => {
                let s = x.sin();
                let c = x.cos();
                c.div(&s)
            }
            Self::Sec => Jet::one(deg).div(&x.cos()),
            Self::Csc => Jet::one(deg).div(&x.sin()),
            Self::Sinh => x.sinh(),
            Self::Cosh => x.cosh(),
            Self::Tanh => {
                let s = x.sinh();
                let c = x.cosh();
                s.div(&c)
            }
            Self::Coth => {
                let s = x.sinh();
                let c = x.cosh();
                c.div(&s)
            }
            Self::Sech => Jet::one(deg).div(&x.cosh()),
            Self::Csch => Jet::one(deg).div(&x.sinh()),

            Self::Sqrt => x.sqrt(),
            Self::Cbrt => {
                // cbrt(x) = exp((1/3)·ln(x))
                let one_third = Cmplx::new(s(1.0 / 3.0), scalar_zero());
                x.ln().scale(&one_third).exp()
            }

            // --- Special: Sinc = sin(x) / x ---
            Self::Sinc => {
                if x.coeffs[0].is_zero() {
                    // sinc(0) = 1 − x²/6 + x⁴/120 − …  (even series)
                    let mut coeffs = vec![Cmplx::zero(); deg + 1];
                    coeffs[0] = Cmplx::one();
                    for n in (2..=deg).step_by(2) {
                        let sign = if ((n >> 1) & 1) == 0 {
                            Cmplx::one()
                        } else {
                            Cmplx::one().neg()
                        };
                        // Taylor coeff: (−1)^{n/2} / (n+1)!
                        // But sinc[n] = f^(n)(0)/n! where f^(n)(0) = (−1)^{n/2}·n!/(n+1)! = (−1)^{n/2}/(n+1)
                        // Actually, sinc(x) = 1 − x²/6 + x⁴/120 − x⁶/5040 + …
                        // coeffs[2] = -1/6  = -1/(3!)     … wait: 3! = 6, so coeffs[2] = −1/6 = −1/3!
                        // coeffs[4] =  1/120 = 1/5!
                        // coeffs[6] = −1/5040 = −1/7!
                        // So coeffs[n] = (−1)^{n/2} / (n+1)! for even n, 0 for odd n
                        let fact = (2..=(n + 1)).fold(Cmplx::one(), |acc, i| {
                            acc.mul(&Cmplx::new(s(i), scalar_zero()))
                        });
                        coeffs[n] = sign.div(&fact);
                    }
                    Jet::new(coeffs)
                } else {
                    x.sin().div(x)
                }
            }

            // --- Var / Const / algebraic ---
            Self::Var => x.clone(),
            Self::Const(c) => Jet::constant(c.clone(), deg),
            Self::Neg(a) => a.eval_jet(x).neg(),
            Self::Inv(a) => a.eval_jet(x).inv(),
            Self::Add(args) => {
                let mut acc = Jet::zero(deg);
                for a in args {
                    acc = acc.add(&a.eval_jet(x));
                }
                acc
            }
            Self::Mul(args) => {
                let mut acc = Jet::one(deg);
                for a in args {
                    acc = acc.mul(&a.eval_jet(x));
                }
                acc
            }
            Self::Compose(outer, inner) => {
                let inner_jet = inner.eval_jet(x);
                outer.eval_jet(&inner_jet)
            }

            // --- Functions with non-circular derivatives ---
            //
            // These functions' derivative ASTs are Jet-compatible (no reference
            // back to the function itself), so the derivative formula works.
            Self::Asin
            | Self::Acos
            | Self::Atan
            | Self::Acot
            | Self::Asec
            | Self::Acsc
            | Self::Asinh
            | Self::Acosh
            | Self::Atanh
            | Self::Acoth
            | Self::Asech
            | Self::Acsch
            | Self::LGamma
            | Self::Erf
            | Self::Erfc => self.noncircular_eval_jet(x),

            // --- EllipticE and EllipticK have circular derivatives ---
            //
            // K'(k) = (E − (1−k²)·K) / (k·(1−k²))
            // E'(k) = (E − K) / k
            //
            // Both derivatives contain E/K, so the derivative formula would
            // recurse infinitely.  Compute both Jets simultaneously via the
            // coupled recurrence.
            Self::EllipticE | Self::EllipticK => {
                let (ej, kj) = Self::elliptic_ek_pair(x);
                if matches!(self, Self::EllipticE) {
                    ej
                } else {
                    kj
                }
            }

            // --- Functions whose derivative is the same family shifted ---
            //
            // Polygamma(n).derivative() = Polygamma(n+1)
            // ZetaDeriv(n).derivative() = ZetaDeriv(n+1)
            // Zeta.derivative() = ZetaDeriv(1)
            //
            // The fallback would recurse infinitely, so we compute the Taylor
            // coefficients directly.  This is valid because these functions
            // only appear as top-level AST nodes (never inside Compose), so
            // `x` is always a simple variable Jet.
            //
            // For a simple variable Jet x = λ + ε:
            //   coeffs[n] = f^(n)(λ) / n!
            Self::Polygamma(m) => Self::polygamma_jet(*m, x),
            Self::Zeta => Self::zeta_jet(x),
            Self::ZetaDeriv(m) => Self::zetaderiv_jet(*m, x),

            // --- Gamma has circular derivative (contains Gamma itself) ---
            //
            // Uses: Γ(x)' = Γ(x) · ψ(x) · x'   (ψ = digamma)
            // For a simple variable Jet:  γ[n] = (1/n)·Σ_{j=0}^{n-1} γ[j]·ψ[n-1-j]
            Self::Gamma => Self::gamma_jet(x),
        }
    }

    /// Direct Jet for Polygamma(m) via Taylor expansion at x₀.
    /// Only correct for simple variable inputs (`Jet::variable`).
    fn polygamma_jet(m: usize, x: &Jet) -> Jet {
        let deg = x.degree();
        let mut coeffs = vec![Cmplx::zero(); deg + 1];
        for (n, c) in coeffs.iter_mut().enumerate() {
            let p = Self::Polygamma(m + n).eval(&x.coeffs[0]);
            let fact = (2..=n)
                .map(|i| Cmplx::new(s(i), scalar_zero()))
                .fold(Cmplx::one(), |acc, f| acc.mul(&f));
            *c = if n == 0 { p } else { p.div(&fact) };
        }
        Jet::new(coeffs)
    }

    /// Direct Jet for Zeta via Taylor expansion at x₀.
    /// Only correct for simple variable inputs (`Jet::variable`).
    fn zeta_jet(x: &Jet) -> Jet {
        let deg = x.degree();
        let mut coeffs = vec![Cmplx::zero(); deg + 1];
        coeffs[0] = Self::Zeta.eval(&x.coeffs[0]);
        for (n, c) in coeffs.iter_mut().enumerate().skip(1) {
            let p = Self::ZetaDeriv(n).eval(&x.coeffs[0]);
            let fact = (2..=n)
                .map(|i| Cmplx::new(s(i), scalar_zero()))
                .fold(Cmplx::one(), |acc, f| acc.mul(&f));
            *c = p.div(&fact);
        }
        Jet::new(coeffs)
    }

    /// Direct Jet for ZetaDeriv(m) via Taylor expansion at x₀.
    /// Only correct for simple variable inputs (`Jet::variable`).
    fn zetaderiv_jet(m: usize, x: &Jet) -> Jet {
        let deg = x.degree();
        let mut coeffs = vec![Cmplx::zero(); deg + 1];
        for (n, c) in coeffs.iter_mut().enumerate() {
            let p = Self::ZetaDeriv(m + n).eval(&x.coeffs[0]);
            let fact = (2..=n)
                .map(|i| Cmplx::new(s(i), scalar_zero()))
                .fold(Cmplx::one(), |acc, f| acc.mul(&f));
            *c = if n == 0 { p } else { p.div(&fact) };
        }
        Jet::new(coeffs)
    }

    /// Direct Jet for Gamma via coupled digamma recurrence.
    /// Only correct for simple variable inputs (`Jet::variable`).
    ///
    /// Uses: Γ(x)' = Γ(x) · ψ(x), so for Jet coefficients:
    ///   γ[0] = Γ(x₀)
    ///   γ[n] = (1/n)·Σ_{j=0}^{n-1} γ[j] · `ψ_jet.coeffs`[n-1-j]
    fn gamma_jet(x: &Jet) -> Jet {
        let deg = x.degree();
        let psi_jet = Self::polygamma_jet(0, x);
        let mut coeffs = vec![Cmplx::zero(); deg + 1];
        coeffs[0] = Self::Gamma.eval(&x.coeffs[0]);
        for n in 1..=deg {
            let mut sum = Cmplx::zero();
            #[allow(
                clippy::needless_range_loop,
                reason = "j indexes coeffs[j] and psi_jet.coeffs[n-1-j] — asymmetric access pattern"
            )]
            for j in 0..n {
                sum = sum.add(&coeffs[j].mul(&psi_jet.coeffs[n - 1 - j]));
            }
            coeffs[n] = sum.div(&Cmplx::new(s(n), scalar_zero()));
        }
        Jet::new(coeffs)
    }

    /// Compute both `EllipticE` and `EllipticK` Jets simultaneously via the
    /// coupled ODE system:
    ///
    ///   E'(k) = (E(k) − K(k)) / k
    ///   K'(k) = (E(k) − (1−k²)·K(k)) / (k·(1−k²))
    ///
    /// Only correct for simple variable inputs (`Jet::variable`).
    fn elliptic_ek_pair(x: &Jet) -> (Jet, Jet) {
        let deg = x.degree();
        let inv_x = x.inv();
        let one_m_x2 = Jet::one(deg).sub(&x.mul(x));
        let scale_k = inv_x.mul(&one_m_x2.inv());

        let mut e = vec![Cmplx::zero(); deg + 1];
        let mut k = vec![Cmplx::zero(); deg + 1];
        e[0] = Self::EllipticE.eval(&x.coeffs[0]);
        k[0] = Self::EllipticK.eval(&x.coeffs[0]);

        // For simple variable Jet x = λ + ε:
        //   e[n+1] = (1/(n+1))·Σ_{j=0}^{n} (e[j]−k[j])·inv_x[n−j]
        //   k[n+1] = (1/(n+1))·Σ_{j=0}^{n} ((e − (1−k²)·k)[j])·scale_k[n−j]
        for n in 0..deg {
            let mut sum_e = Cmplx::zero();
            let mut sum_k = Cmplx::zero();
            for j in 0..=n {
                let diff = e[j].sub(&k[j]);
                sum_e = sum_e.add(&diff.mul(&inv_x.coeffs[n - j]));

                let mut k_num = e[j].clone();
                for i in 0..=j {
                    k_num = k_num.sub(&one_m_x2.coeffs[i].mul(&k[j - i]));
                }
                sum_k = sum_k.add(&k_num.mul(&scale_k.coeffs[n - j]));
            }
            let np1 = Cmplx::new(s(n + 1), scalar_zero());
            e[n + 1] = sum_e.div(&np1);
            k[n + 1] = sum_k.div(&np1);
        }
        (Jet::new(e), Jet::new(k))
    }
    /// Fallback for functions that lack a direct Jet recurrence but whose
    /// derivative AST IS Jet-compatible (no circular reference back to self).
    /// Uses derivative relation: f'(x) = df(x), then recovers f's Taylor
    /// coefficients from df's and x's via the shift identity:
    ///   (n)·coeffs[n] = Σ_{i=0}^{n-1} `df_jet`[i] · (n-i) · x.coeffs[n-i]
    fn noncircular_eval_jet(&self, x: &Jet) -> Jet {
        let deg = x.degree();
        let mut coeffs = vec![Cmplx::zero(); deg + 1];
        coeffs[0] = self.eval(&x.coeffs[0]);

        if deg == 0 {
            return Jet::new(coeffs);
        }

        let df = self.derivative();
        // If derivative is NaN (Erf, EllipticK/E), no higher coefficients
        if matches!(&df, Self::Const(c) if !c.0.is_finite() || !c.1.is_finite()) {
            return Jet::new(coeffs);
        }

        let df_jet = df.eval_jet(x);
        for (n, c) in coeffs.iter_mut().enumerate().skip(1) {
            let mut sum = Cmplx::zero();
            #[allow(
                clippy::needless_range_loop,
                reason = "i indexes df_jet.coeffs[i] and x.coeffs[n-i]"
            )]
            for i in 0..n {
                let weight = Cmplx::new(s(n - i), scalar_zero());
                sum = sum.add(&df_jet.coeffs[i].mul(&weight).mul(&x.coeffs[n - i]));
            }
            *c = sum.div(&Cmplx::new(s(n), scalar_zero()));
        }
        Jet::new(coeffs)
    }
}

// ============================================================================
// Helper: apply an AST function via spectral decomposition
// ============================================================================

fn apply_via_spectral(mv: &CliffordNumber, ast: &SpectralFn) -> CliffordNumber {
    let gens = mv.generator_set();

    if is_pure_scalar(mv) {
        let s = ast.eval(&Cmplx::new(mv.coeff(0).clone(), scalar_zero())).0;
        let mut out = CliffordNumber::zero_unchecked(gens.clone());
        out.coeffs_mut_slice()[0] = s;
        return out;
    }

    // Path 1: ≤3 generators via Mat(2,ℂ) spectral decomposition
    if let Some(blade_mats) = compute_blade_basis(gens) {
        let mat = to_matrix_with_basis(mv, &blade_mats);

        // Standard Eigendecomposition
        if let Some((u, u_inv, l1, l2)) = eigendecompose(&mat) {
            let fd = Mat2C::new(ast.eval(&l1), Cmplx::zero(), Cmplx::zero(), ast.eval(&l2));
            let result_mat = u.mul(&fd).mul(&u_inv);
            return from_matrix_with_gens(&result_mat, gens, &blade_mats);
        }

        // Defective matrix fallback: compute via Schur form T = Q^H * A * Q
        let (q_mat, t_mat) = schur_decompose_2x2(&mat);
        let lambda = t_mat.a.clone();
        let f_lambda = ast.eval(&lambda);
        let f_prime_lambda = ast.derivative().eval(&lambda);
        if f_lambda.0.is_finite() && f_prime_lambda.0.is_finite() {
            let fd = Mat2C::new(
                f_lambda.clone(),
                t_mat.b.mul(&f_prime_lambda),
                Cmplx::zero(),
                f_lambda,
            );
            let q_dagger = q_mat.dagger();
            let result_mat = q_mat.mul(&fd).mul(&q_dagger);
            return from_matrix_with_gens(&result_mat, gens, &blade_mats);
        }
        return nan_clifford(gens);
    }

    // Path 2: Fallback to large_mat.rs for large or highly nilpotent algebras
    large_apply_via_spectral(mv, ast)
}

fn extract_scalar(mv: &CliffordNumber) -> Scalar {
    let blade0 = mv.coeff(0);
    if blade0.is_zero() {
        scalar_zero()
    } else {
        blade0.clone()
    }
}

fn is_pure_scalar(mv: &CliffordNumber) -> bool {
    for blade in 1..mv.blade_count() {
        if !mv.coeff(blade).is_zero() {
            return false;
        }
    }
    true
}

/// Apply a real-scalar binary function `f(order, x)` to a `CliffordNumber`.
fn apply_real_binary_fn<F>(mv: &CliffordNumber, order: &CliffordNumber, f: &F) -> CliffordNumber
where
    F: Fn(&Scalar, &Scalar) -> Scalar,
{
    if !is_pure_scalar(order) {
        return nan_clifford(&mv.gens);
    }
    let n = extract_scalar(order);
    if is_pure_scalar(mv) {
        let s = f(mv.coeff(0), &n);
        let mut out = CliffordNumber::zero_unchecked(mv.gens.clone());
        out.coeffs_mut_slice()[0] = s;
        out
    } else if can_embed(&mv.gens) {
        apply_via_spectral_closure(mv, &|c: &Cmplx| c.apply_real_only(&|s| f(s, &n)))
    } else {
        nan_clifford(&mv.gens)
    }
}

// Old fallback for binary functions (since they don't have AST nodes yet)
fn apply_via_spectral_closure<F>(mv: &CliffordNumber, f: &F) -> CliffordNumber
where
    F: Fn(&Cmplx) -> Cmplx,
{
    let gens = mv.generator_set();
    if let Some(blade_mats) = compute_blade_basis(gens) {
        let mat = to_matrix_with_basis(mv, &blade_mats);
        if let Some((u, u_inv, l1, l2)) = eigendecompose(&mat) {
            let fd = Mat2C::new(f(&l1), Cmplx::zero(), Cmplx::zero(), f(&l2));
            let result_mat = u.mul(&fd).mul(&u_inv);
            return from_matrix_with_gens(&result_mat, gens, &blade_mats);
        }
        // Defective matrix fallback via Schur form + numeric derivative
        let (q_mat, t_mat) = schur_decompose_2x2(&mat);
        let lambda = t_mat.a.clone();
        let f_lambda = f(&lambda);
        if f_lambda.0.is_finite() && f_lambda.1.is_zero() {
            // Numeric derivative via central finite differences
            let eps = Scalar::epsilon();
            let h_cmplx = Cmplx::new(eps.clone(), scalar_zero());
            let f_plus = f(&lambda.add(&h_cmplx));
            let f_minus = f(&lambda.sub(&h_cmplx));
            let two_h = Cmplx::new(s(2.0) * eps, scalar_zero());
            let f_prime_lambda = (f_plus.sub(&f_minus)).div(&two_h);
            if f_prime_lambda.0.is_finite() && f_prime_lambda.1.is_zero() {
                let fd = Mat2C::new(
                    f_lambda.clone(),
                    t_mat.b.mul(&f_prime_lambda),
                    Cmplx::zero(),
                    f_lambda,
                );
                let q_dagger = q_mat.dagger();
                let result_mat = q_mat.mul(&fd).mul(&q_dagger);
                return from_matrix_with_gens(&result_mat, gens, &blade_mats);
            }
        }
        return nan_clifford(gens);
    }
    // Fallback to large_mat.rs
    large_apply_via_spectral_closure(mv, f)
}

pub fn nan_clifford(gens: &GeneratorSet) -> CliffordNumber {
    let mut mv = CliffordNumber::zero_unchecked(gens.clone());
    mv.coeffs_mut_slice()[0] = scalar_nan();
    mv
}

/// Compute the Schur decomposition of a 2×2 matrix.
/// Returns `(Q, T)` where `Q` is unitary and `T = Q^H * A * Q` is upper triangular.
/// For a defective matrix, T = [[λ, x], [0, λ]] with repeated eigenvalue λ.
fn schur_decompose_2x2(mat: &Mat2C) -> (Mat2C, Mat2C) {
    // Eigenvalue: λ = trace / 2 (for 2x2, repeated eigenvalue case)
    let trace = mat.a.add(&mat.d);
    let half = Cmplx::new(s(0.5), scalar_zero());
    let lambda = trace.mul(&half);

    // Compute eigenvector v = [v1, v2] for λ
    let v1 = Cmplx::one();
    let v2 = if mat.b.is_zero() && mat.c.is_zero() {
        Cmplx::zero()
    } else if !mat.b.is_zero() {
        let diff = lambda.sub(&mat.a);
        diff.div(&mat.b)
    } else {
        // c·v₁ + (d−λ)·v₂ = 0, with v₁ = 1 → v₂ = c / (λ−d)
        let diff = lambda.sub(&mat.d);
        mat.c.div(&diff)
    };

    // Normalize v
    let v_norm = v1.abs_sq().add(&v2.abs_sq()).sqrt();
    let v_norm_c = Cmplx::new(v_norm, scalar_zero());
    let v1n = v1.div(&v_norm_c);
    let v2n = v2.div(&v_norm_c);

    // Find orthogonal vector w = [-conj(v2), conj(v1)]
    let w1 = v2n.conj().neg();
    let w2 = v1n.conj();

    // Build Q = [v, w]
    let q_mat = Mat2C::new(v1n, w1, v2n, w2);

    // Compute T = Q^H * A * Q
    let qh = q_mat.dagger();
    let qh_a = qh.mul(mat);
    let t_mat = qh_a.mul(&q_mat);

    (q_mat, t_mat)
}

// ============================================================================
// Geometric inverse via spectral decomposition
// ============================================================================

impl CliffordNumber {
    /// Computes the geometric inverse of the Clifford number.
    #[must_use]
    pub fn geometric_inverse(&self) -> Self {
        let gens = self.generator_set();

        // Path 1: ≤3 generators via Mat(2,ℂ)
        if let Some(blade_mats) = compute_blade_basis(gens) {
            let mat = to_matrix_with_basis(self, &blade_mats);
            if let Some((u, u_inv, l1, l2)) = eigendecompose(&mat) {
                let one = Cmplx::one();
                let inv_l1 = one.div(&l1);
                let inv_l2 = one.div(&l2);
                let inv_d = Mat2C::new(inv_l1, Cmplx::zero(), Cmplx::zero(), inv_l2);
                let result_mat = u.mul(&inv_d).mul(&u_inv);
                return from_matrix_with_gens(&result_mat, gens, &blade_mats);
            }
            return nan_clifford(gens);
        }

        // Path 2: Fallback to large_mat.rs
        apply_via_spectral(self, &SpectralFn::Inv(Box::new(SpectralFn::Var)))
    }
}

// ============================================================================
// Number trait implementation for CliffordNumber
// ============================================================================

impl Number for CliffordNumber {
    fn sin(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Sin)
    }
    fn cos(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Cos)
    }
    fn tan(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Tan)
    }
    fn cot(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Cot)
    }
    fn sec(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Sec)
    }
    fn csc(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Csc)
    }

    fn asin(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Asin)
    }
    fn acos(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Acos)
    }
    fn atan(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Atan)
    }
    fn acot(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Acot)
    }
    fn asec(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Asec)
    }
    fn acsc(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Acsc)
    }

    fn sinh(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Sinh)
    }
    fn cosh(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Cosh)
    }
    fn tanh(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Tanh)
    }
    fn coth(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Coth)
    }
    fn sech(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Sech)
    }
    fn csch(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Csch)
    }

    fn asinh(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Asinh)
    }
    fn acosh(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Acosh)
    }
    fn atanh(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Atanh)
    }
    fn acoth(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Acoth)
    }
    fn asech(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Asech)
    }
    fn acsch(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Acsch)
    }

    fn exp(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Exp)
    }
    fn expm1(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Expm1)
    }
    fn exp_neg(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::ExpNeg)
    }
    fn ln(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Ln)
    }
    fn log1p(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Log1p)
    }

    fn sqrt(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Sqrt)
    }
    fn cbrt(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Cbrt)
    }

    fn erf(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Erf)
    }
    fn erfc(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Erfc)
    }
    fn gamma(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Gamma)
    }
    fn lgamma(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::LGamma)
    }
    fn digamma(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Polygamma(0))
    }
    fn trigamma(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Polygamma(1))
    }
    fn tetragamma(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Polygamma(2))
    }

    fn sinc(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Sinc)
    }
    fn elliptic_k(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::EllipticK)
    }
    fn elliptic_e(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::EllipticE)
    }
    fn zeta(&self) -> Self {
        apply_via_spectral(self, &SpectralFn::Zeta)
    }
    fn exp_polar(&self) -> Self {
        self.exp()
    }

    fn atan2(&self, x: &Self) -> Self {
        if is_pure_scalar(self) && is_pure_scalar(x) {
            let s = self.coeff(0).atan2(x.coeff(0));
            let mut mv = Self::zero_unchecked(self.gens.clone());
            mv.coeffs_mut_slice()[0] = s;
            mv
        } else {
            nan_clifford(&self.gens)
        }
    }

    fn log_base(&self, base: &Self) -> Self {
        if is_pure_scalar(self) && is_pure_scalar(base) {
            let s = self.coeff(0).log_base(base.coeff(0));
            let mut mv = Self::zero_unchecked(self.gens.clone());
            mv.coeffs_mut_slice()[0] = s;
            mv
        } else {
            nan_clifford(&self.gens)
        }
    }

    fn pow(&self, exp: &Self) -> Self {
        if is_pure_scalar(self) && is_pure_scalar(exp) {
            let s = self.coeff(0).pow(exp.coeff(0));
            let mut mv = Self::zero_unchecked(self.gens.clone());
            mv.coeffs_mut_slice()[0] = s;
            mv
        } else if is_pure_scalar(exp) && can_embed(&self.gens) {
            let exp_s = extract_scalar(exp);
            let f = |c: &Cmplx| {
                let log_c = c.ln();
                log_c.scale(&exp_s).exp()
            };
            large_apply_via_spectral_closure(self, &f)
        } else {
            nan_clifford(&self.gens)
        }
    }

    fn besselj(&self, order: &Self) -> Self {
        apply_real_binary_fn(self, order, &|s, n| s.besselj(n))
    }
    fn bessely(&self, order: &Self) -> Self {
        apply_real_binary_fn(self, order, &|s, n| s.bessely(n))
    }
    fn besseli(&self, order: &Self) -> Self {
        apply_real_binary_fn(self, order, &|s, n| s.besseli(n))
    }
    fn besselk(&self, order: &Self) -> Self {
        apply_real_binary_fn(self, order, &|s, n| s.besselk(n))
    }
    fn polygamma(&self, order: &Self) -> Self {
        apply_real_binary_fn(self, order, &|s, n| s.polygamma(n))
    }

    fn beta(&self, other: &Self) -> Self {
        if is_pure_scalar(self) && is_pure_scalar(other) {
            let s = self.coeff(0).beta(other.coeff(0));
            let mut mv = Self::zero_unchecked(self.gens.clone());
            mv.coeffs_mut_slice()[0] = s;
            mv
        } else {
            nan_clifford(&self.gens)
        }
    }

    fn zeta_deriv(&self, order: &Self) -> Self {
        apply_real_binary_fn(self, order, &|s, n| s.zeta_deriv(n))
    }

    fn lambertw(&self, n: &Self) -> Self {
        apply_real_binary_fn(self, n, &|s, o| s.lambertw(o))
    }
    fn hermite(&self, n: &Self) -> Self {
        apply_real_binary_fn(self, n, &|s, o| s.hermite(o))
    }
    fn assoc_legendre(&self, _l: &Self, _m: &Self) -> Self {
        nan_clifford(&self.gens)
    }
    fn spherical_harmonic(&self, _l: &Self, _m: &Self, _phi: &Self) -> Self {
        nan_clifford(&self.gens)
    }
    fn is_zero(&self) -> bool {
        for i in 0..self.blade_count() {
            if !self.coeff(i).is_zero() {
                return false;
            }
        }
        true
    }
    fn is_one(&self) -> bool {
        if !self.coeff(0).is_one() {
            return false;
        }
        for i in 1..self.blade_count() {
            if !self.coeff(i).is_zero() {
                return false;
            }
        }
        true
    }
    fn is_neg_one(&self) -> bool {
        if !self.coeff(0).is_neg_one() {
            return false;
        }
        for i in 1..self.blade_count() {
            if !self.coeff(i).is_zero() {
                return false;
            }
        }
        true
    }
    fn is_integer(&self) -> bool {
        is_pure_scalar(self) && self.coeff(0).is_integer()
    }
    fn is_negative(&self) -> bool {
        is_pure_scalar(self) && self.coeff(0).is_negative()
    }
    fn is_positive(&self) -> bool {
        is_pure_scalar(self) && self.coeff(0).is_positive()
    }
    fn is_finite(&self) -> bool {
        for i in 0..self.blade_count() {
            if !self.coeff(i).is_finite() {
                return false;
            }
        }
        true
    }
    fn to_float(&self) -> Self {
        let mut out = self.clone();
        for i in 0..out.blade_count() {
            out.coeffs_mut_slice()[i] = out.coeff(i).to_float();
        }
        out
    }
    fn approx_eq_number(&self, other: &Self, tolerance: &Self) -> bool {
        if self.gens != other.gens {
            return false;
        }
        let tol = extract_scalar(tolerance);
        for i in 0..self.blade_count() {
            if !self.coeff(i).approx_eq_number(other.coeff(i), &tol) {
                return false;
            }
        }
        true
    }
    fn total_cmp(&self, other: &Self) -> Ordering {
        if self.gens != other.gens {
            return self.gens.len().cmp(&other.gens.len());
        }
        let limit = self.blade_count();
        for i in 0..limit {
            let cmp = self.coeff(i).total_cmp(other.coeff(i));
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        Ordering::Equal
    }
    fn num_max(&self, other: &Self) -> Self {
        if self.total_cmp(other) == Ordering::Less {
            other.clone()
        } else {
            self.clone()
        }
    }
    fn num_min(&self, other: &Self) -> Self {
        if self.total_cmp(other) == Ordering::Greater {
            other.clone()
        } else {
            self.clone()
        }
    }

    fn abs(&self) -> Self {
        if is_pure_scalar(self) {
            let mut out = self.clone();
            out.coeffs_mut_slice()[0] = self.coeff(0).abs();
            out
        } else {
            nan_clifford(&self.gens)
        }
    }
    fn signum(&self) -> Self {
        if is_pure_scalar(self) {
            let mut out = self.clone();
            out.coeffs_mut_slice()[0] = self.coeff(0).signum();
            out
        } else {
            nan_clifford(&self.gens)
        }
    }
    fn floor(&self) -> Self {
        if is_pure_scalar(self) {
            let mut out = self.clone();
            out.coeffs_mut_slice()[0] = self.coeff(0).floor();
            out
        } else {
            nan_clifford(&self.gens)
        }
    }
    fn ceil(&self) -> Self {
        if is_pure_scalar(self) {
            let mut out = self.clone();
            out.coeffs_mut_slice()[0] = self.coeff(0).ceil();
            out
        } else {
            nan_clifford(&self.gens)
        }
    }
    fn round(&self) -> Self {
        if is_pure_scalar(self) {
            let mut out = self.clone();
            out.coeffs_mut_slice()[0] = self.coeff(0).round();
            out
        } else {
            nan_clifford(&self.gens)
        }
    }
    fn fract(&self) -> Self {
        if is_pure_scalar(self) {
            let mut out = self.clone();
            out.coeffs_mut_slice()[0] = self.coeff(0).fract();
            out
        } else {
            nan_clifford(&self.gens)
        }
    }
    fn negate(&self) -> Self {
        -self
    }
}
