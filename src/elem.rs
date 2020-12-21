// Copyright 2020 Yao Pengfei.
//
// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
// SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
// OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
// CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

use crate::curve::{
    bn_add_mod, bn_mont_pro, bn_point_add, bn_point_mul, bn_scalar_add_mod, bn_scalar_mod,
    bn_scalar_mont_pro, bn_scalar_sub_mod, bn_scalar_to_inv, bn_to_inv,
};
use crate::param::CurveCtx;
use crate::public::{Point, PublicKey, BN_LENGTH};
use core::marker::PhantomData;
use num_bigint::BigUint;
use num_traits::{One, Zero};

// Indicates that the element is not encoded; there is no *R* factor
// that needs to be canceled out.
#[derive(Copy, Clone)]
pub enum Unencoded {}

// Indicates that the element is encoded; the value has one *R*
// factor that needs to be canceled out.
#[derive(Copy, Clone)]
pub enum R {}

// Indicates the element is encoded twice; the value has two *R*
// factors that need to be canceled out.
#[derive(Copy, Clone)]
pub enum RR {}

// Indicates the element is inversely encoded; the value has one
// 1/*R* factor that needs to be canceled out.
#[derive(Copy, Clone)]
pub enum RInverse {}

pub trait Encoding {}

impl Encoding for RR {}
impl Encoding for R {}
impl Encoding for Unencoded {}
impl Encoding for RInverse {}

/// The encoding of the result of a reduction.
pub trait ReductionEncoding {
    type Output: Encoding;
}

impl ReductionEncoding for RR {
    type Output = R;
}
impl ReductionEncoding for R {
    type Output = Unencoded;
}
impl ReductionEncoding for Unencoded {
    type Output = RInverse;
}

/// The encoding of the result of a multiplication.
pub trait ProductEncoding {
    type Output: Encoding;
}

impl<E: ReductionEncoding> ProductEncoding for (Unencoded, E) {
    type Output = E::Output;
}

impl<E: Encoding> ProductEncoding for (R, E) {
    type Output = E;
}

impl<E: ReductionEncoding> ProductEncoding for (RInverse, E)
where
    E::Output: ReductionEncoding,
{
    type Output = <<E as ReductionEncoding>::Output as ReductionEncoding>::Output;
}

// XXX: Rust doesn't allow overlapping impls,
// TODO (if/when Rust allows it):
// impl<E1, E2: ReductionEncoding> ProductEncoding for
//         (E1, E2) {
//     type Output = <(E2, E1) as ProductEncoding>::Output;
// }
impl ProductEncoding for (RR, Unencoded) {
    type Output = <(Unencoded, RR) as ProductEncoding>::Output;
}
impl ProductEncoding for (RR, RInverse) {
    type Output = <(RInverse, RR) as ProductEncoding>::Output;
}

/// Elements are always fully reduced with respect to *m*; i.e.
/// the 0 <= x < m for every value x.
#[derive(Clone)]
pub struct Elem<M> {
    pub inner: BigUint,

    /// The modulus *m* for the ring ℤ/mℤ for which this element is a value.
    pub m: PhantomData<M>,
}

impl<M> Elem<M> {
    pub fn zero() -> Self {
        Self {
            inner: BigUint::zero(),
            m: PhantomData,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.inner == BigUint::zero()
    }

    pub fn is_equal(&self, other: &Elem<M>) -> bool {
        self.inner == other.inner
    }

    pub fn bytes_less_safe(&self) -> [u8; BN_LENGTH] {
        let mut r = [0; BN_LENGTH];
        r.copy_from_slice(&self.inner.to_bytes_be());
        r
    }
}

pub fn elem_mul<EA: Encoding, EB: Encoding>(
    a: &Elem<EA>,
    b: &Elem<EB>,
    cctx: &CurveCtx,
) -> Elem<<(EA, EB) as ProductEncoding>::Output>
where
    (EA, EB): ProductEncoding,
{
    Elem {
        inner: bn_mont_pro(&a.inner, &b.inner, cctx),
        m: PhantomData,
    }
}

pub fn elem_add(a: &Elem<R>, b: &Elem<R>, cctx: &CurveCtx) -> Elem<R> {
    Elem {
        inner: bn_add_mod(&a.inner, &b.inner, cctx),
        m: PhantomData,
    }
}

pub fn elem_inv_sqr_to_mont(a: &Elem<R>, cctx: &CurveCtx) -> Elem<R> {
    assert!(!a.inner.is_zero());
    let a_inv = bn_to_inv(&a.inner, cctx);

    Elem {
        inner: bn_mont_pro(&a_inv, &a_inv, cctx),
        m: PhantomData,
    }
}

pub fn elem_to_unencoded(a: &Elem<R>, cctx: &CurveCtx) -> Elem<Unencoded> {
    Elem {
        inner: bn_mont_pro(&a.inner, &BigUint::one(), cctx),
        m: PhantomData,
    }
}

pub fn elem_reduced_to_scalar(e: &Elem<Unencoded>, cctx: &CurveCtx) -> Scalar {
    Scalar {
        inner: bn_scalar_mod(&e.inner, cctx),
        m: PhantomData,
    }
}

pub fn scalar_to_elem(e: &Scalar) -> Elem<Unencoded> {
    Elem {
        inner: e.inner.clone(),
        m: PhantomData,
    }
}

/// A scalar. Its value is in [0, n). Zero-valued scalars are forbidden in most
/// contexts.
pub type Scalar<N = Unencoded> = Elem<N>;

pub fn scalar_inv_to_mont(a: &Scalar, cctx: &CurveCtx) -> Scalar<R> {
    Scalar {
        inner: bn_scalar_to_inv(&a.inner, cctx),
        m: PhantomData,
    }
}

pub fn scalar_to_unencoded(a: &Scalar<R>, cctx: &CurveCtx) -> Scalar {
    Scalar {
        inner: bn_scalar_mont_pro(&a.inner, &BigUint::one(), cctx),
        m: PhantomData,
    }
}

pub fn scalar_mul<EA: Encoding, EB: Encoding>(
    a: &Scalar<EA>,
    b: &Scalar<EB>,
    cctx: &CurveCtx,
) -> Scalar<<(EA, EB) as ProductEncoding>::Output>
where
    (EA, EB): ProductEncoding,
{
    Scalar {
        inner: bn_scalar_mont_pro(&a.inner, &b.inner, cctx),
        m: PhantomData,
    }
}

pub fn scalar_add(a: &Scalar, b: &Scalar, cctx: &CurveCtx) -> Scalar {
    Scalar {
        inner: bn_scalar_add_mod(&a.inner, &b.inner, cctx),
        m: PhantomData,
    }
}

pub fn scalar_sub(a: &Scalar, b: &Scalar, cctx: &CurveCtx) -> Scalar {
    Scalar {
        inner: bn_scalar_sub_mod(&a.inner, &b.inner, cctx),
        m: PhantomData,
    }
}

pub fn scalar_g(g_scalar: &Scalar, cctx: &CurveCtx) -> [BigUint; 3] {
    bn_point_mul(&cctx.g_point, &g_scalar.inner, cctx)
}

pub fn scalar_p(p_scalar: &Scalar, pk: &PublicKey, cctx: &CurveCtx) -> [BigUint; 3] {
    bn_point_mul(&pk.to_point(cctx).to_bns(), &p_scalar.inner, cctx)
}

pub fn twin_mul(g_scalar: &Scalar, p_scalar: &Scalar, pk: &PublicKey, cctx: &CurveCtx) -> Point {
    let g_point = scalar_g(g_scalar, cctx);
    let p_point = scalar_p(p_scalar, pk, cctx);
    Point::new(bn_point_add(&g_point, &p_point, cctx))
}
