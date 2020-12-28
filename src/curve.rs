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

use crate::param::CurveCtx;
use num_bigint::BigUint;
use num_traits::identities::Zero;
use num_traits::One;

const CURVE_LENGTH: usize = 256;

pub(crate) fn bn_add_mod(a: &BigUint, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH && b.bits() <= CURVE_LENGTH);
    (a + b) % &cctx.p
}

fn bn_neg_mod(a: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH);
    let a = a % &cctx.p;
    &cctx.p - &a
}

fn bn_sub_mod(a: &BigUint, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH && b.bits() <= CURVE_LENGTH);
    let neg_b = bn_neg_mod(b, cctx);
    (a + &neg_b) % &cctx.p
}

pub(crate) fn bn_mul_mod(a: &BigUint, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH);
    bn_mod(&(b * a), cctx)
}

// a << b
fn bn_shl_mod(a: &BigUint, b: usize, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH);
    let a_shl = a << b;
    bn_mod(&a_shl, cctx)
}

fn bn_mod(a: &BigUint, cctx: &CurveCtx) -> BigUint {
    let mut rem = a.clone();
    let mut rem_len = rem.bits();
    while rem_len > CURVE_LENGTH {
        let left_rem = &rem >> CURVE_LENGTH;
        rem -= (&left_rem + (&left_rem >> 32) + (&left_rem >> 160) - (&left_rem >> 192)) * &cctx.p;
        rem_len = rem.bits();
    }
    if rem >= cctx.p {
        rem -= &cctx.p;
    }
    rem
}

/// The algorithm: "add-1998-cmo-2"
/// Cost: 12M + 4S + 6add + 1*2.
///       Z1Z1 = Z1^2
///       Z2Z2 = Z2^2
///       U1 = X1*Z2Z2
///       U2 = X2*Z1Z1
///       S1 = Y1*Z2*Z2Z2
///       S2 = Y2*Z1*Z1Z1
///       H = U2-U1
///       HH = H^2
///       HHH = H*H^H
///       r = S2-S1
///       V = U1*HH
///       X3 = r^2-HHH-2*V
///       Y3 = r*(V-X3)-S1*HHH
///       Z3 = Z1*Z2*H
pub(crate) fn bn_point_add(a: &[BigUint; 3], b: &[BigUint; 3], cctx: &CurveCtx) -> [BigUint; 3] {
    let a_x = &a[0];
    let a_y = &a[1];
    let a_z = &a[2];
    assert!(a_x.bits() <= CURVE_LENGTH && a_y.bits() <= CURVE_LENGTH && a_z.bits() <= CURVE_LENGTH);
    let b_x = &b[0];
    let b_y = &b[1];
    let b_z = &b[2];
    assert!(b_x.bits() <= CURVE_LENGTH && b_y.bits() <= CURVE_LENGTH && b_z.bits() <= CURVE_LENGTH);

    if a_z.is_zero() {
        return b.clone();
    } else if b_z.is_zero() {
        return a.clone();
    } else if a_x == b_x && a_y == b_y && a_z == b_z {
        let rem_arr = bn_point_double(a, cctx);
        return rem_arr;
    }

    let a_z_sqr = bn_mul_mod(&a_z, &a_z, cctx);
    let b_z_sqr = bn_mul_mod(&b_z, &b_z, cctx);
    let u1 = bn_mul_mod(&a_x, &b_z_sqr, cctx);
    let u2 = bn_mul_mod(&b_x, &a_z_sqr, cctx);
    let a_z_cub = bn_mul_mod(&a_z_sqr, &a_z, cctx);
    let b_z_cub = bn_mul_mod(&b_z_sqr, &b_z, cctx);
    let s1 = bn_mul_mod(&a_y, &b_z_cub, cctx);
    let s2 = bn_mul_mod(&b_y, &a_z_cub, cctx);
    let h = bn_sub_mod(&u2, &u1, cctx);
    let r = bn_sub_mod(&s2, &s1, cctx);
    let r_sqr = bn_mul_mod(&r, &r, cctx);
    let h_sqr = bn_mul_mod(&h, &h, cctx);
    let h_cub = bn_mul_mod(&h_sqr, &h, cctx);

    let lam1 = bn_mul_mod(&u1, &h_sqr, cctx); // u1*h^2
    let rem_x = bn_sub_mod(
        &bn_sub_mod(&r_sqr, &h_cub, cctx),
        &bn_shl_mod(&lam1, 1, cctx),
        cctx,
    );
    let rem_y = bn_sub_mod(
        &bn_mul_mod(&r, &bn_sub_mod(&lam1, &rem_x, cctx), cctx),
        &bn_mul_mod(&s1, &h_cub, cctx),
        cctx,
    );
    let rem_z = bn_mul_mod(&bn_mul_mod(&a_z, &b_z, cctx), &h, cctx);

    [rem_x, rem_y, rem_z]
}

/// The algorithm: "dbl-2001-b"
/// Cost: 3M + 5S + 8add + 1*3 + 1*4 + 2*8
///      delta = Z12
///       gamma = Y12
///       beta = X1*gamma
///       alpha = 3*(X1-delta)*(X1+delta)
///       X3 = alpha2-8*beta
///       Z3 = (Y1+Z1)2-gamma-delta
///       Y3 = alpha*(4*beta-X3)-8*gamma2
pub(crate) fn bn_point_double(a: &[BigUint; 3], cctx: &CurveCtx) -> [BigUint; 3] {
    let a_x = &a[0];
    let a_y = &a[1];
    let a_z = &a[2];
    assert!(a_x.bits() <= CURVE_LENGTH && a_y.bits() <= CURVE_LENGTH && a_z.bits() <= CURVE_LENGTH);
    let delta = bn_mul_mod(a_z, a_z, cctx);
    let gamma = bn_mul_mod(a_y, a_y, cctx);
    let beta = bn_mul_mod(a_x, &gamma, cctx);
    let alpha = bn_mul_mod(
        &bn_mul_mod(
            &bn_sub_mod(a_x, &delta, cctx),
            &bn_add_mod(a_x, &delta, cctx),
            cctx,
        ),
        &BigUint::from(3 as usize),
        cctx,
    );
    let rem_x = bn_sub_mod(
        &bn_mul_mod(&alpha, &alpha, cctx),
        &bn_shl_mod(&beta, 3, cctx),
        cctx,
    );
    let lam1 = bn_sub_mod(&bn_shl_mod(&beta, 2, cctx), &rem_x, cctx); // 4 * beta - x3
    let rem_y = bn_sub_mod(
        &bn_mul_mod(&alpha, &lam1, cctx),
        &bn_shl_mod(&bn_mul_mod(&gamma, &gamma, cctx), 3, cctx),
        cctx,
    );
    let lam2 = bn_add_mod(a_y, a_z, cctx);
    let rem_z = bn_sub_mod(
        &bn_sub_mod(&bn_mul_mod(&lam2, &lam2, cctx), &gamma, cctx),
        &delta,
        cctx,
    );
    [rem_x, rem_y, rem_z]
}

pub(crate) fn bn_point_mul(a: &[BigUint; 3], scalar: &BigUint, cctx: &CurveCtx) -> [BigUint; 3] {
    let a_x = &a[0];
    let a_y = &a[1];
    let a_z = &a[2];
    assert!(
        a_x.bits() <= CURVE_LENGTH
            && a_y.bits() <= CURVE_LENGTH
            && a_z.bits() <= CURVE_LENGTH
            && scalar.bits() <= CURVE_LENGTH
    );
    let scalar_bz = scalar.to_bytes_le();
    let mut a_order = [a[0].clone(), a[1].clone(), a[2].clone()];
    let mut rem = [BigUint::zero(), BigUint::zero(), BigUint::zero()];

    for scalar_byte in scalar_bz {
        let mut bit: usize = 0;
        while bit < 8 {
            if (scalar_byte >> bit) & 0x01 != 0 {
                rem = bn_point_add(&rem, &a_order, cctx);
            }
            a_order = bn_point_double(&a_order, cctx);
            bit += 1;
        }
    }
    rem
}

// (`a` squared `squarings` times) * b
#[inline]
pub(crate) fn bn_sqr_mul(a: &BigUint, squarings: usize, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(squarings >= 1 && a.bits() <= CURVE_LENGTH);
    let mut rem = bn_mul_mod(a, a, cctx);
    for _ in 1..squarings {
        rem = bn_mul_mod(&rem, &rem, cctx);
    }
    bn_mul_mod(&rem, b, cctx)
}

pub(crate) fn bn_to_inv(a: &BigUint, cctx: &CurveCtx) -> BigUint {
    // Calculate the modular inverse of scalar |a| using Fermat's Little
    // Theorem:
    // Calculate a**-1 (mod q) == a**(q - 2) (mod q)
    //
    // The exponent (p - 2) is:
    //
    //    0xfffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffd

    assert!(a.bits() <= CURVE_LENGTH);
    let b_1 = a;
    let b_11 = bn_sqr_mul(b_1, 1, b_1, cctx);
    let b_111 = bn_sqr_mul(&b_11, 1, b_1, cctx);
    let f_11 = bn_sqr_mul(&b_111, 3, &b_111, cctx);
    let fff = bn_sqr_mul(&f_11, 6, &f_11, cctx);
    let fff_111 = bn_sqr_mul(&fff, 3, &b_111, cctx);
    let fffffff_11 = bn_sqr_mul(&fff_111, 15, &fff_111, cctx);
    let ffffffff = bn_sqr_mul(&fffffff_11, 2, &b_11, cctx);

    // fffffff_111
    let mut acc = bn_sqr_mul(&fffffff_11, 1, &b_1, cctx);

    // fffffffe
    acc = bn_mul_mod(&acc, &acc, cctx);

    // fffffffeffffffff
    acc = bn_sqr_mul(&acc, 32, &ffffffff, cctx);

    // fffffffeffffffffffffffff
    acc = bn_sqr_mul(&acc, 32, &ffffffff, cctx);

    // fffffffeffffffffffffffffffffffff
    acc = bn_sqr_mul(&acc, 32, &ffffffff, cctx);

    // fffffffeffffffffffffffffffffffffffffffff
    acc = bn_sqr_mul(&acc, 32, &ffffffff, cctx);

    // fffffffeffffffffffffffffffffffffffffffff00000000ffffffff
    acc = bn_sqr_mul(&acc, 64, &ffffffff, cctx);

    // fffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffff_11
    acc = bn_sqr_mul(&acc, 30, &fffffff_11, cctx);

    // fffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffd
    bn_sqr_mul(&acc, 2, b_1, cctx)
}

#[cfg(test)]
fn bn_to_affine(a: &[BigUint; 3], cctx: &CurveCtx) -> [BigUint; 2] {
    let a_x = &a[0];
    let a_y = &a[1];
    let a_z = &a[2];
    assert!(a_x.bits() <= CURVE_LENGTH && a_y.bits() <= CURVE_LENGTH && a_z.bits() <= CURVE_LENGTH);

    let z_inv = bn_to_inv(&a_z, cctx);
    let zz_inv = bn_mul_mod(&z_inv, &z_inv, cctx);
    let zzz_inv = bn_mul_mod(&zz_inv, &z_inv, cctx);

    let rem_x = bn_mul_mod(&a_x, &zz_inv, cctx);
    let rem_y = bn_mul_mod(&a_y, &zzz_inv, cctx);

    [rem_x, rem_y]
}

pub(crate) fn bn_to_jacobi(a: &[BigUint; 2]) -> [BigUint; 3] {
    // 1 * r modsm2p256
    [a[0].clone(), a[1].clone(), BigUint::one()]
}

pub(crate) fn bn_scalar_mul_mod(a: &BigUint, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH);
    a * b % &cctx.n
}

pub(crate) fn bn_scalar_add_mod(a: &BigUint, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH && b.bits() <= CURVE_LENGTH);
    (a + b) % &cctx.n
}

fn bn_scalar_neg_mod(a: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH);
    let a = a % &cctx.n;
    &cctx.n - &a
}

pub(crate) fn bn_scalar_sub_mod(a: &BigUint, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(a.bits() <= CURVE_LENGTH && b.bits() <= CURVE_LENGTH);
    let neg_b = bn_scalar_neg_mod(b, cctx);
    (a + &neg_b) % &cctx.n
}

// `a` squared `squarings` times * `b`
fn bn_scalar_sqr_mul(a: &BigUint, squarings: usize, b: &BigUint, cctx: &CurveCtx) -> BigUint {
    assert!(squarings >= 1 && a.bits() <= CURVE_LENGTH);
    let mut rem = bn_scalar_mul_mod(a, a, cctx);
    for _ in 1..squarings {
        rem = bn_scalar_mul_mod(&rem, &rem, cctx);
    }
    bn_scalar_mul_mod(&rem, b, cctx)
}

pub(crate) fn bn_scalar_to_inv(a: &BigUint, cctx: &CurveCtx) -> BigUint {
    // Calculate the modular inverse of scalar |a| using Fermat's Little
    // Theorem:
    //
    //    a**-1 (mod n) == a**(n - 2) (mod n)
    //
    // The exponent (n - 2) is:
    //
    //    0xfffffffeffffffffffffffffffffffff7203df6b21c6052b53bbf40939d54121

    // Indexes into `d`.
    const B_1: usize = 0;
    const B_10: usize = 1;
    const B_11: usize = 2;
    const B_101: usize = 3;
    const B_111: usize = 4;
    const B_1111: usize = 5;
    const B_10101: usize = 6;
    const B_101111: usize = 7;

    let mut d = [
        BigUint::zero(),
        BigUint::zero(),
        BigUint::zero(),
        BigUint::zero(),
        BigUint::zero(),
        BigUint::zero(),
        BigUint::zero(),
        BigUint::zero(),
    ];

    d[B_1] = a.clone();
    d[B_10] = bn_scalar_mul_mod(&d[B_1], &d[B_1], cctx);
    d[B_11] = bn_scalar_mul_mod(&d[B_10], &d[B_1], cctx);
    d[B_101] = bn_scalar_mul_mod(&d[B_10], &d[B_11], cctx);
    d[B_111] = bn_scalar_mul_mod(&d[B_101], &d[B_10], cctx);
    let b_1010 = bn_scalar_mul_mod(&d[B_101], &d[B_101], cctx);
    d[B_1111] = bn_scalar_mul_mod(&b_1010, &d[B_101], cctx);
    d[B_10101] = bn_scalar_sqr_mul(&b_1010, 0 + 1, &d[B_1], cctx);
    let b_101010 = bn_scalar_mul_mod(&d[B_10101], &d[B_10101], cctx);
    d[B_101111] = bn_scalar_mul_mod(&b_101010, &d[B_101], cctx);
    let b_111111 = bn_scalar_mul_mod(&b_101010, &d[B_10101], cctx);
    let b_1111111 = bn_scalar_sqr_mul(&b_111111, 0 + 1, &d[B_1], cctx);

    let ff = bn_scalar_sqr_mul(&b_111111, 0 + 2, &d[B_11], cctx);
    let ffff = bn_scalar_sqr_mul(&ff, 0 + 8, &ff, cctx);
    let ffffffff = bn_scalar_sqr_mul(&ffff, 0 + 16, &ffff, cctx);

    // ffffff
    let mut acc = bn_scalar_sqr_mul(&ffff, 0 + 8, &ff, cctx);

    // fffffff_111
    acc = bn_scalar_sqr_mul(&acc, 0 + 7, &b_1111111, cctx);

    // fffffffe
    acc = bn_scalar_mul_mod(&acc, &acc, cctx);

    // fffffffeffffffff
    acc = bn_scalar_sqr_mul(&acc, 0 + 32, &ffffffff, cctx);

    // fffffffeffffffffffffffff
    acc = bn_scalar_sqr_mul(&acc, 0 + 32, &ffffffff, cctx);

    // fffffffeffffffffffffffffffffffff
    acc = bn_scalar_sqr_mul(&acc, 0 + 32, &ffffffff, cctx);

    // The rest of the exponent, in binary, is:
    //
    //    0111,001,00000001111,01111,101,10101,1,001,0000111,00011,000000101,0010101,1
    //    111,1,00111,0111,00111,0010101,1,0000101111,11,00011,00011,001,0010101,001111

    //    0111,001,00000001111,01111,101,10101,1,001,0000111,00011,000000101,0010101,
    //    10101,00111,0111,01111,11,01,0000001,001,00111,00111,010101,01,000001,001,00001

    static REMAINING_WINDOWS: [(usize, usize); 27] = [
        (1 + 3, B_111),
        (2 + 1, B_1),
        (7 + 4, B_1111),
        (1 + 4, B_1111),
        (0 + 3, B_101),
        (0 + 5, B_10101),
        (0 + 1, B_1),
        (2 + 1, B_1),
        (4 + 3, B_111),
        (3 + 2, B_11),
        (6 + 3, B_101),
        (2 + 5, B_10101),
        (0 + 5, B_10101),
        (2 + 3, B_111),
        (1 + 3, B_111),
        (1 + 4, B_1111),
        (0 + 2, B_11),
        (1 + 1, B_1),
        (6 + 1, B_1),
        (2 + 1, B_1),
        (2 + 3, B_111),
        (2 + 3, B_111),
        (1 + 5, B_10101),
        (1 + 1, B_1),
        (5 + 1, B_1),
        (2 + 1, B_1),
        (4 + 1, B_1),
    ];

    for &(squarings, digit) in &REMAINING_WINDOWS {
        acc = bn_scalar_sqr_mul(&acc, squarings, &d[digit], cctx);
    }

    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    #[test]
    fn bn_mul_mod_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                .unwrap(),
        );
        println!("bn_mul_mod_test 1: {:x}", &bn_mul_mod(&a, &a, cctx));

        println!(
            "bn_mul_mod_test 2: a * 1: {:x}, a: {:x}",
            bn_mul_mod(&a, &cctx.r_p, cctx),
            a
        );
    }

    #[test]
    fn bn_neg_mod_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffffeffffffffffffffffffffffffffffffff00000001ffffffffffffffff")
                .unwrap(),
        );
        println!("bn_neg_mod_test: {:x}", bn_neg_mod(&a, cctx));
    }

    #[test]
    fn bn_sub_mod_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                .unwrap(),
        );
        let b = BigUint::from_bytes_be(
            &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                .unwrap(),
        );
        println!("bn_sub_mod_test 1: {:x}", bn_sub_mod(&a, &b, cctx));

        let a = BigUint::from_bytes_be(
            &hex::decode("0100000000000000000000000000000001000000000000000000000001").unwrap(),
        );
        let b = BigUint::from_bytes_be(
            &hex::decode("fffffffeffffffffffffffffffffffffffffffff00000001ffffffffffffffff")
                .unwrap(),
        );
        println!("bn_sub_mod_test 2: {:x}", bn_sub_mod(&a, &b, cctx));
    }

    #[test]
    fn bn_to_inv_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                .unwrap(),
        );
        let res_mon = bn_to_inv(&a, cctx);
        println!(
            "bn_to_inv_test: {:x}",
            bn_mul_mod(&res_mon, &BigUint::new(vec![1]), cctx)
        );
    }

    #[test]
    fn bn_shl_mod_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                .unwrap(),
        );
        println!("bn_shl_mod_test: {:x}", bn_shl_mod(&a, 7, cctx));
    }

    #[test]
    fn bn_point_double_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let ori_point_g = [
            BigUint::from_bytes_be(
                &hex::decode("32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7")
                    .unwrap(),
            ),
            BigUint::from_bytes_be(
                &hex::decode("bc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0")
                    .unwrap(),
            ),
        ];
        let projective_mont_point_g = bn_to_jacobi(&ori_point_g);
        let double_projective_mont_point_g = bn_point_double(&projective_mont_point_g, cctx);
        println!(
            "bn_point_double_test 1: x: {:x}, y: {:x}",
            double_projective_mont_point_g[0], double_projective_mont_point_g[1]
        );
        let double_affine_mont_point_g = bn_to_affine(&double_projective_mont_point_g, cctx);
        println!(
            "bn_point_double_test 2: double g_x: {:x}, g_y: {:x}",
            double_affine_mont_point_g[0], double_affine_mont_point_g[1]
        );
    }

    #[test]
    fn bn_point_add_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let g_2 = [
            BigUint::from_bytes_be(
                &hex::decode("0d7e9c18caa5736a5349d94b5788cd2483bdc9ba2d8fa9380af037bfbc3be46a")
                    .unwrap(),
            ),
            BigUint::from_bytes_be(
                &hex::decode("947e74656c21bdf5c7b145169b7157acccbd8d37c4a8e82b6a7e1a1d69db9ac1")
                    .unwrap(),
            ),
        ];
        let g_4 = [
            BigUint::from_bytes_be(
                &hex::decode("50dc8e3ac899dbe18a86bcb4a09f9020487ea27fe9016209393f7c5a98615060")
                    .unwrap(),
            ),
            BigUint::from_bytes_be(
                &hex::decode("6ffc31c525bce9e34d0bd55632cf70ed1de135ea7c7383bdfc099043fd619998")
                    .unwrap(),
            ),
        ];
        let pro_g_2 = bn_to_jacobi(&g_2);
        let pro_g_4 = bn_to_jacobi(&g_4);
        let pro_g_6 = bn_point_add(&pro_g_2, &pro_g_4, cctx);
        println!(
            "bn_point_add_test 1: x: {:x}, y: {:x}",
            pro_g_6[0], pro_g_6[1]
        );
        let aff_g_6 = bn_to_affine(&pro_g_6, cctx);
        println!(
            "bn_point_add_test 2: 6g_x: {:x}, 6g_y: {:x}",
            aff_g_6[0], aff_g_6[1]
        )
    }

    #[test]
    fn bn_point_mul_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let ori_point_g = [
            BigUint::from_bytes_be(
                &hex::decode("32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7")
                    .unwrap(),
            ),
            BigUint::from_bytes_be(
                &hex::decode("bc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0")
                    .unwrap(),
            ),
        ];
        let projective_mont_point_g = bn_to_jacobi(&ori_point_g);

        let scalar = bn_shl_mod(&BigUint::new(vec![31]), 8 * 1, cctx);
        let pro_point = bn_point_mul(&projective_mont_point_g, &scalar, cctx);
        println!(
            "bn_point_mul_test 1: x: {:x}, y: {:x}",
            pro_point[0], pro_point[1]
        );
        let aff_point = bn_to_affine(&pro_point, cctx);
        println!(
            "bn_point_mul_test 2: affine_point: {:x}, affine_point: {:x}",
            aff_point[0], aff_point[1]
        )
    }

    #[test]
    fn bn_scalar_mont_pro_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(&hex::decode("01").unwrap());
        let b = BigUint::from_bytes_be(
            &hex::decode("fffffffeffffffffffffffffffffffff7203df6b21c6052b53bbf40939d54122")
                .unwrap(),
        );
        println!(
            "bn_scalar_mont_pro_test: {:x}",
            bn_scalar_mul_mod(&a, &b, cctx)
        );
    }

    #[test]
    fn bn_scalar_sqr_mul_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                .unwrap(),
        );
        let b = BigUint::from_bytes_be(
            &hex::decode("010000000000000000000000008dfc2094de39fad4ac440bf6c62abedd").unwrap(),
        );
        println!(
            "bn_scalar_sqr_mul_test: {:x}",
            bn_scalar_sqr_mul(&a, 5, &b, cctx)
        );
    }

    #[test]
    fn bn_sqr_mul_test() {
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                .unwrap(),
        );
        let b = BigUint::from_bytes_be(
            &hex::decode("52ab139ac09ec8307bdb6926ab664658d3f55c3f46cdfd7516553623adc0a99a")
                .unwrap(),
        );
        println!("bn_sqr_mul_test: {:x}", bn_sqr_mul(&a, 4, &b, cctx));
    }
}

#[cfg(feature = "internal_benches")]
mod internal_benches {
    use super::*;
    use num_bigint::BigUint;
    extern crate test;

    #[bench]
    fn bn_to_inv_bench(bench: &mut test::Bencher) {
        // This benchmark assumes that `elem_inverse_squared()` is
        // constant-time so inverting 1 mod q is as good of a choice as
        // anything.
        let cctx = &CurveCtx::sm2p256_new();
        let a = BigUint::from_bytes_be(&hex::decode("01").unwrap());
        bench.iter(|| {
            let _ = bn_to_inv(&a, cctx);
        });
    }
}
