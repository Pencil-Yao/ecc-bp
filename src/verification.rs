use crate::curve::{bn_mul_mod, bn_point_add, bn_point_mul, bn_scalar_add_mod, bn_scalar_sub_mod};
use crate::error::KeyRejected;
use crate::exchange::verify_jacobian_point_is_on_the_curve;
use crate::param::CurveCtx;
use crate::public::{Point, PublicKey, BN_LENGTH};
use num_bigint::BigUint;

pub struct Signature {
    r: BigUint,
    s: BigUint,
}

impl Signature {
    pub fn new(r: &[u8; BN_LENGTH], s: &[u8; BN_LENGTH]) -> Result<Self, KeyRejected> {
        let r = BigUint::from_bytes_be(r);
        let s = BigUint::from_bytes_be(s);

        Ok(Signature { r, s })
    }

    pub fn from_scalars(r: BigUint, s: BigUint) -> Self {
        Signature { r, s }
    }

    pub fn r(&self) -> [u8; BN_LENGTH] {
        let mut r_out = [0; BN_LENGTH];
        r_out.copy_from_slice(&self.r.to_bytes_be()[..BN_LENGTH]);
        r_out
    }

    pub fn s(&self) -> [u8; BN_LENGTH] {
        let mut s_out = [0; BN_LENGTH];
        s_out.copy_from_slice(&self.s.to_bytes_be()[..BN_LENGTH]);
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
        let e = BigUint::from_bytes_be(digest) % &cctx.n;

        let u = bn_scalar_add_mod(&self.r, &self.s, cctx);
        let r = bn_scalar_sub_mod(&self.r, &e, cctx);

        let g_point = bn_point_mul(&cctx.g_point, &self.s, cctx);
        let p_point = bn_point_mul(&pk.to_point().to_bns(), &u, cctx);
        let point = Point::new(bn_point_add(&g_point, &p_point, cctx));

        verify_jacobian_point_is_on_the_curve(&point, cctx)?;

        fn sig_r_equals_x(r: &BigUint, point: &Point, cctx: &CurveCtx) -> bool {
            let x = point.point_x();
            let z = point.point_z();
            let z2 = bn_mul_mod(&z, &z, cctx);
            let r_jacobian = bn_mul_mod(&z2, &r, cctx);
            r_jacobian == x
        }

        if sig_r_equals_x(&r, &point, cctx) {
            return Ok(());
        }
        Err(KeyRejected::verify_digest_error())
    }
}
