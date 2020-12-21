use crate::elem::{
    elem_mul, elem_reduced_to_scalar, elem_to_unencoded, scalar_add, scalar_sub, scalar_to_elem,
    twin_mul, Elem, Scalar, Unencoded,
};
use crate::error::KeyRejected;
use crate::exchange::verify_jacobian_point_is_on_the_curve;
use crate::param::CurveCtx;
use crate::public::{Point, PublicKey, BN_LENGTH};
use num_bigint::BigUint;
use std::marker::PhantomData;

pub struct Signature {
    r: Scalar,
    s: Scalar,
}

impl Signature {
    pub fn new(r: &[u8; BN_LENGTH], s: &[u8; BN_LENGTH]) -> Result<Self, KeyRejected> {
        let rl = BigUint::from_bytes_be(r);
        let r = Scalar {
            inner: rl,
            m: PhantomData,
        };

        let sl = BigUint::from_bytes_be(s);
        let s = Scalar {
            inner: sl,
            m: PhantomData,
        };

        Ok(Signature { r, s })
    }

    pub fn from_scalars(r: Scalar, s: Scalar) -> Self {
        Signature { r, s }
    }

    pub fn r(&self) -> [u8; BN_LENGTH] {
        let mut r_out = [0; BN_LENGTH];
        r_out.copy_from_slice(&self.r.bytes_less_safe());
        r_out
    }

    pub fn s(&self) -> [u8; BN_LENGTH] {
        let mut s_out = [0; BN_LENGTH];
        s_out.copy_from_slice(&self.s.bytes_less_safe());
        s_out
    }

    pub fn sm2_verify(
        &self,
        pk: &PublicKey,
        msg: &[u8],
        cctx: &CurveCtx,
    ) -> Result<(), KeyRejected> {
        let digest = (cctx.hasher)(&pk, msg)?;

        self.sm2_verify_digest(pk, &digest, cctx)
    }

    pub fn sm2_verify_digest(
        &self,
        pk: &PublicKey,
        digest: &[u8],
        cctx: &CurveCtx,
    ) -> Result<(), KeyRejected> {
        let dl = BigUint::from_bytes_be(digest);
        let edl = Elem {
            inner: dl,
            m: PhantomData,
        };
        let e = elem_reduced_to_scalar(&edl, cctx);

        let (u1, u2) = (&self.s, scalar_add(&self.r, &self.s, cctx));
        let r = scalar_sub(&self.r, &e, cctx);

        let point = twin_mul(&u1, &u2, &pk, cctx);

        verify_jacobian_point_is_on_the_curve(&point, cctx)?;

        fn sig_r_equals_x(r: &Elem<Unencoded>, point: &Point, cctx: &CurveCtx) -> bool {
            let x = point.point_x();
            let z = point.point_z();
            let z2 = elem_mul(&z, &z, cctx);
            let r_jacobian = elem_mul(&z2, &r, cctx);
            let x = elem_to_unencoded(&x, cctx);
            r_jacobian.is_equal(&x)
        }

        let r = scalar_to_elem(&r);
        if sig_r_equals_x(&r, &point, cctx) {
            return Ok(());
        }
        Err(KeyRejected::verify_digest_error())
    }
}
