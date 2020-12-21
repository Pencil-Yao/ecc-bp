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
use crate::public::BN_LENGTH;
use crate::rand::SecureRandom;
use crate::signing::KeyPair;
use crate::verification::Signature;
use ::rand::{prelude::ThreadRng, thread_rng, Rng};

mod curve;
mod elem;
mod error;
mod exchange;
mod param;
mod private;
mod public;
mod rand;
mod signing;
mod verification;

fn main() {
    pub struct EgRand(ThreadRng);

    impl SecureRandom for EgRand {
        fn fill(&mut self, dest: &mut [u8]) {
            self.0.fill(dest)
        }
    }

    let test_word = b"hello world";
    let mut rng = EgRand(thread_rng());
    let cctx = CurveCtx::sm2p256_new();

    let mut private_key = [0; BN_LENGTH];
    rng.fill(&mut private_key);

    private_key.copy_from_slice(
        &hex::decode("b8aa2a5bd9a9cf448984a247e63cb3878859d02b886e1bc63cd5c6dd46a744ab").unwrap(),
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
