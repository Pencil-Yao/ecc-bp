use crate::curve::{bn_add_mod, bn_mul_mod, bn_to_inv};
use crate::error::KeyRejected;
use crate::param::CurveCtx;
use crate::public::Point;
use num_bigint::BigUint;

pub fn big_endian_affine_from_jacobian(
    point: &Point,
    cctx: &CurveCtx,
) -> Result<(Vec<u8>, Vec<u8>), KeyRejected> {
    let (x, y) = affine_from_jacobian(&point, cctx)?;

    Ok((x.to_bytes_be(), y.to_bytes_be()))
}

pub fn affine_from_jacobian(
    point: &Point,
    cctx: &CurveCtx,
) -> Result<(BigUint, BigUint), KeyRejected> {
    let x = point.point_x();
    let y = point.point_y();
    let z = point.point_z();

    let zz_inv = bn_to_inv(&z, cctx);

    let x_aff = bn_mul_mod(&x, &zz_inv, cctx);

    let y_aff = {
        let zzzz_inv = bn_mul_mod(&zz_inv, &zz_inv, cctx);
        let zzz_inv = bn_mul_mod(&z, &zzzz_inv, cctx);
        bn_mul_mod(&y, &zzz_inv, cctx)
    };

    verify_affine_point_is_on_the_curve((&x_aff, &y_aff), &cctx.a, &cctx.b, cctx)?;

    Ok((x_aff, y_aff))
}

pub fn verify_jacobian_point_is_on_the_curve(
    point: &Point,
    cctx: &CurveCtx,
) -> Result<(), KeyRejected> {
    let z = point.point_z();

    if z == BigUint::from(0 as usize) {
        return Err(KeyRejected::zero_error());
    }

    let x = point.point_x();
    let y = point.point_y();

    let z2 = bn_mul_mod(&z, &z, cctx);
    let z4 = bn_mul_mod(&z2, &z2, cctx);
    let z4_a = bn_mul_mod(&z4, &cctx.a, cctx);
    let z6 = bn_mul_mod(&z4, &z2, cctx);
    let z6_b = bn_mul_mod(&z6, &cctx.b, cctx);

    verify_affine_point_is_on_the_curve((&x, &y), &z4_a, &z6_b, cctx)
}

pub fn verify_affine_point_is_on_the_curve(
    (x, y): (&BigUint, &BigUint),
    a: &BigUint,
    b: &BigUint,
    cctx: &CurveCtx,
) -> Result<(), KeyRejected> {
    let lhs = bn_mul_mod(y, y, cctx);

    let x2 = bn_mul_mod(x, x, cctx);
    let x2_a = bn_add_mod(&x2, a, cctx);
    let x2_a_x = bn_mul_mod(&x2_a, x, cctx);
    let rhs = bn_add_mod(&x2_a_x, b, cctx);

    if lhs != rhs {
        return Err(KeyRejected::not_on_curve_error());
    }
    Ok(())
}
