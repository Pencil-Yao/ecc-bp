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
    bn_point_mul, bn_scalar_add_mod, bn_scalar_mul_mod, bn_scalar_sub_mod, bn_scalar_to_inv,
};
use crate::error::KeyRejected;
use crate::exchange::{affine_from_jacobian, big_endian_affine_from_jacobian};
use crate::param::CurveCtx;
use crate::private::create_private_key;
use crate::public::{Point, PublicKey, BN_LENGTH};
use crate::rand::SecureRandom;
use crate::verification::Signature;
use num_bigint::BigUint;
use num_traits::{One, Zero};

pub struct KeyPair {
    d: BigUint,
    pk: PublicKey,
}

impl KeyPair {
    pub fn new(private_key: &[u8; BN_LENGTH], cctx: &CurveCtx) -> Result<Self, KeyRejected> {
        let d = BigUint::from_bytes_be(private_key) % &cctx.n;
        let pk = public_from_private(&d, cctx)?;
        Ok(KeyPair { d, pk })
    }

    pub fn public_key(&self) -> PublicKey {
        self.pk
    }

    pub fn sm2_sign(
        &self,
        rng: &mut dyn SecureRandom,
        message: &[u8],
        cctx: &CurveCtx,
    ) -> Result<Signature, KeyRejected> {
        let digest = (cctx.hasher)(&self.public_key(), &message)?;

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
            #[cfg(test)]
            let rk = BigUint::from_bytes_be(
                &hex::decode("fffffc4d0000064efffffb8c00000324fffffdc600000543fffff8950000053b")
                    .unwrap(),
            );

            let rq = Point::new(bn_point_mul(&cctx.g_point, &rk, cctx));

            let r = {
                let (x, _) = affine_from_jacobian(&rq, cctx)?;
                x % &cctx.n
            };
            if r.is_zero() {
                continue;
            }

            let e = BigUint::from_bytes_be(&digest);

            let r = bn_scalar_add_mod(&r, &e, cctx);

            let left = bn_scalar_to_inv(&bn_scalar_add_mod(&self.d, &BigUint::one(), cctx), cctx);
            let dr = bn_scalar_mul_mod(&self.d, &r, cctx);
            let right = bn_scalar_sub_mod(&rk, &dr, cctx);
            let s = bn_scalar_mul_mod(&left, &right, cctx);

            #[cfg(test)]
            {
                let target_r = BigUint::from_bytes_be(
                    &hex::decode("80511be00b753e05b0b7abe51be3753b17151244aa5f66e6f87939d3e00a3b4d")
                        .unwrap(),
                );
                let target_s = BigUint::from_bytes_be(
                    &hex::decode("39286762f8a6cd0ab02b66d99c9a76cc49bf92c3b0f64aff9d80eaf350bb37bc")
                        .unwrap(),
                );
                assert_eq!(r, target_r);
                assert_eq!(s, target_s)
            }

            return Ok(Signature::from_scalars(r, s));
        }
        Err(KeyRejected::sign_digest_error())
    }
}

pub fn public_from_private(
    private_key: &BigUint,
    cctx: &CurveCtx,
) -> Result<PublicKey, KeyRejected> {
    let pk_point = Point::new(bn_point_mul(&cctx.g_point, private_key, cctx));

    let (x, y) = big_endian_affine_from_jacobian(&pk_point, cctx)?;

    Ok(PublicKey::new(&x, &y))
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

        sig2.sm2_verify(&key_pair.public_key(), test_word, &cctx)
            .unwrap()
    }
}
