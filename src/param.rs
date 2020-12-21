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

use crate::elem::{Elem, R};
use crate::error::KeyRejected;
use crate::public::{PublicKey, BN_LENGTH};
use num_bigint::BigUint;
use std::marker::PhantomData;

/// ecc equation: y^2 == x^3 +ax + b (modp)
/// r = 2 ^256
pub struct CurveCtx {
    // r -1
    pub r_sub_one: BigUint,
    pub p: BigUint,
    // Calculate the modular inverse of scalar |a| using Fermat's Little
    // Theorem:
    //
    //    a**-1 (mod n) == a**(n - 2) (mod n)
    // so, p_inv_r_neg == -p^-1 modr, ps. p^-1 = p^(r-2) modr
    pub p_inv_r_neg: BigUint,
    // r modp
    pub r_p: BigUint,
    // r^2 modp
    pub rr_p: BigUint,
    pub n: BigUint,
    // like p_inv_r_neg, n_inv_r_neg == -n^-1 modr, ps. n^-1 = n^(r-2) modr
    pub n_inv_r_neg: BigUint,
    // r^2 modn
    pub rr_n: BigUint,
    // a * r modp
    pub a_mont: Elem<R>,
    // b * r modp
    pub b_mont: Elem<R>,
    // generator point at jacobi
    // x_aff = to_mont(x), y_aff = to_mont(y)
    // g_point = to_jacobi(x_aff, y_aff)
    pub g_point: [BigUint; 3],
    // hash function
    pub hasher: fn(pk: &PublicKey, msg: &[u8]) -> Result<[u8; BN_LENGTH], KeyRejected>,
}

impl CurveCtx {
    pub fn sm2p256_new() -> CurveCtx {
        let r = &BigUint::from_bytes_be(
            &hex::decode("010000000000000000000000000000000000000000000000000000000000000000")
                .unwrap(),
        );
        let p = &BigUint::from_bytes_be(
            &hex::decode("fffffffeffffffffffffffffffffffffffffffff00000000ffffffffffffffff")
                .unwrap(),
        );
        let n = &BigUint::from_bytes_be(
            &hex::decode("fffffffeffffffffffffffffffffffff7203df6b21c6052b53bbf40939d54123")
                .unwrap(),
        );
        let a = BigUint::from_bytes_be(
            &hex::decode("fffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffc")
                .unwrap(),
        );
        let b = BigUint::from_bytes_be(
            &hex::decode("28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93")
                .unwrap(),
        );
        let g_x = BigUint::from_bytes_be(
            &hex::decode("32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7")
                .unwrap(),
        );
        let g_y = BigUint::from_bytes_be(
            &hex::decode("bc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0")
                .unwrap(),
        );

        fn sm2hash(pk: &PublicKey, msg: &[u8]) -> Result<[u8; BN_LENGTH], KeyRejected> {
            let ctx = libsm::sm2::signature::SigCtx::new();
            let pk_point = ctx
                .load_pubkey(pk.bytes_less_safe())
                .map_err(|_| KeyRejected::hash_error())?;
            Ok(ctx.hash("1234567812345678", &pk_point, msg))
        }

        let a_mont: Elem<R> = Elem {
            inner: a * r % p,
            m: PhantomData,
        };
        let b_mont: Elem<R> = Elem {
            inner: b * r % p,
            m: PhantomData,
        };

        let ctx = CurveCtx {
            r_sub_one: r - 1 as u32,
            p: p.clone(),
            p_inv_r_neg: BigUint::from_bytes_be(
                &hex::decode("fffffffc00000001fffffffe00000000ffffffff000000010000000000000001")
                    .unwrap(),
            ),
            r_p: r % p,
            rr_p: r * r % p,
            n: n.clone(),
            n_inv_r_neg: BigUint::from_bytes_be(
                &hex::decode("6f39132f82e4c7bc2b0068d3b08941d4df1e8d34fc8319a5327f9e8872350975")
                    .unwrap(),
            ),
            rr_n: r * r % n,
            a_mont,
            b_mont,
            g_point: [g_x * r % p, g_y * r % p, r % p],
            hasher: sm2hash,
        };

        ctx
    }
}
