use crate::elem::{elem_add, elem_inv_sqr_to_mont, elem_mul, elem_to_unencoded, Elem, R};
use crate::error::KeyRejected;
use crate::param::CurveCtx;
use crate::public::{Point, BN_LENGTH};

pub fn big_endian_affine_from_jacobian(
    point: &Point,
    cctx: &CurveCtx,
) -> Result<([u8; BN_LENGTH], [u8; BN_LENGTH]), KeyRejected> {
    let (x_aff, y_aff) = affine_from_jacobian(&point, cctx)?;
    let x = elem_to_unencoded(&x_aff, cctx);
    let y = elem_to_unencoded(&y_aff, cctx);

    Ok((x.bytes_less_safe(), y.bytes_less_safe()))
}

pub fn affine_from_jacobian(
    point: &Point,
    cctx: &CurveCtx,
) -> Result<(Elem<R>, Elem<R>), KeyRejected> {
    let x = point.point_x();
    let y = point.point_y();
    let z = point.point_z();

    let zz_inv = elem_inv_sqr_to_mont(&z, cctx);

    let x_aff = elem_mul(&x, &zz_inv, cctx);

    let y_aff = {
        let zzzz_inv = elem_mul(&zz_inv, &zz_inv, cctx);
        let zzz_inv = elem_mul(&z, &zzzz_inv, cctx);
        elem_mul(&y, &zzz_inv, cctx)
    };

    verify_affine_point_is_on_the_curve((&x_aff, &y_aff), &cctx.a_mont, &cctx.b_mont, cctx)?;

    Ok((x_aff, y_aff))
}

pub fn verify_jacobian_point_is_on_the_curve(
    point: &Point,
    cctx: &CurveCtx,
) -> Result<(), KeyRejected> {
    let z = point.point_z();

    if z.is_zero() {
        return Err(KeyRejected::zero_error());
    }

    let x = point.point_x();
    let y = point.point_y();

    let z2 = elem_mul(&z, &z, cctx);
    let z4 = elem_mul(&z2, &z2, cctx);
    let z4_a = elem_mul(&z4, &cctx.a_mont, cctx);
    let z6 = elem_mul(&z4, &z2, cctx);
    let z6_b = elem_mul(&z6, &cctx.b_mont, cctx);

    verify_affine_point_is_on_the_curve((&x, &y), &z4_a, &z6_b, cctx)
}

pub fn verify_affine_point_is_on_the_curve(
    (x, y): (&Elem<R>, &Elem<R>),
    a: &Elem<R>,
    b: &Elem<R>,
    cctx: &CurveCtx,
) -> Result<(), KeyRejected> {
    let lhs = elem_mul(y, y, cctx);

    let x2 = elem_mul(x, x, cctx);
    let x2_a = elem_add(&x2, a, cctx);
    let x2_a_x = elem_mul(&x2_a, x, cctx);
    let rhs = elem_add(&x2_a_x, b, cctx);

    if !lhs.is_equal(&rhs) {
        return Err(KeyRejected::not_on_curve_error());
    }
    Ok(())
}
