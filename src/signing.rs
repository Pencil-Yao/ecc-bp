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

use crate::curve::bn_scalar_to_mont;
use crate::elem::{
    elem_reduced_to_scalar, elem_to_unencoded, scalar_add, scalar_g, scalar_inv_to_mont,
    scalar_mul, scalar_sub, scalar_to_unencoded, Elem, Scalar, R,
};
use crate::error::KeyRejected;
use crate::exchange::{affine_from_jacobian, big_endian_affine_from_jacobian};
use crate::param::CurveCtx;
use crate::private::create_private_key;
use crate::public::{Point, PublicKey, BN_LENGTH};
use crate::rand::SecureRandom;
use crate::verification::Signature;
use num_bigint::BigUint;
use num_traits::One;
use std::marker::PhantomData;

pub struct KeyPair {
    d: Scalar<R>, // *R*
}

impl KeyPair {
    pub fn new(private_key: &[u8; BN_LENGTH], cctx: &CurveCtx) -> Result<Self, KeyRejected> {
        let key_bn = BigUint::from_bytes_be(private_key);
        let d = Scalar {
            inner: bn_scalar_to_mont(&key_bn, cctx),
            m: PhantomData,
        };
        Ok(KeyPair { d })
    }

    pub fn public_from_private(&self, cctx: &CurveCtx) -> Result<PublicKey, KeyRejected> {
        let du = scalar_to_unencoded(&self.d, cctx);
        let pk_point = Point::new(scalar_g(&du, cctx));

        let (x, y) = big_endian_affine_from_jacobian(&pk_point, cctx)?;

        Ok(PublicKey::new(&x, &y))
    }

    pub fn sm2_sign(
        &self,
        rng: &mut dyn SecureRandom,
        message: &[u8],
        cctx: &CurveCtx,
    ) -> Result<Signature, KeyRejected> {
        let digest = (cctx.hasher)(&self.public_from_private(cctx)?, &message)?;

        self.sm2_sign_digest(rng, &digest, cctx)
    }

    fn sm2_sign_digest(
        &self,
        rng: &mut dyn SecureRandom,
        digest: &[u8],
        cctx: &CurveCtx,
    ) -> Result<Signature, KeyRejected> {
        for _ in 0..100 {
            let rk = create_private_key(rng, cctx)?;

            let rq = Point::new(scalar_g(&rk, cctx));

            let r = {
                let (x, _) = affine_from_jacobian(&rq, cctx)?;
                let x = elem_to_unencoded(&x, cctx);
                elem_reduced_to_scalar(&x, cctx)
            };
            if r.is_zero() {
                continue;
            }

            let dl = BigUint::from_bytes_be(&digest);
            let edl = Elem {
                inner: dl,
                m: PhantomData,
            };
            let e = elem_reduced_to_scalar(&edl, cctx);

            let scalar_one: Scalar = Scalar {
                inner: BigUint::one(),
                m: PhantomData,
            };

            let r = scalar_add(&r, &e, cctx);

            let da_ue = scalar_to_unencoded(&self.d, cctx);
            let left = scalar_inv_to_mont(&scalar_add(&da_ue, &scalar_one, cctx), cctx);
            let dr = scalar_mul(&self.d, &r, cctx);
            let right = scalar_sub(&rk, &dr, cctx);
            let s = scalar_mul(&left, &right, cctx);

            return Ok(Signature::from_scalars(r, s));
        }
        Err(KeyRejected::sign_digest_error())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::ThreadRng;
    use rand::Rng;

    #[test]
    fn sm2_sign_verify_test() {
        pub struct EgRand(ThreadRng);

        impl SecureRandom for EgRand {
            fn fill(&mut self, dest: &mut [u8]) {
                self.0.fill(dest)
            }
        }

        let test_word = b"hello world";
        let mut rng = EgRand(rand::thread_rng());
        let cctx = CurveCtx::sm2p256_new();

        let mut private_key = [0; BN_LENGTH];
        rng.fill(&mut private_key);

        private_key.copy_from_slice(
            &hex::decode("b8aa2a5bd9a9cf448984a247e63cb3878859d02b886e1bc63cd5c6dd46a744ab")
                .unwrap(),
        );

        let key_pair = KeyPair::new(&private_key, &cctx).unwrap();

        let sig = key_pair.sm2_sign(&mut rng, test_word, &cctx).unwrap();

        let r = sig.r();
        let s = sig.s();
        let sig2 = Signature::new(&r, &s).unwrap();

        sig2.sm2_verify(
            &key_pair.public_from_private(&cctx).unwrap(),
            test_word,
            &cctx,
        )
        .unwrap()
    }
}
