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

use crate::curve::{bn_to_jacobi, bn_to_mont};
use crate::elem::{Elem, R};
use crate::param::CurveCtx;
use num_bigint::BigUint;

pub struct Point([BigUint; 3]);

impl Point {
    pub fn new(a: [BigUint; 3]) -> Self {
        Point(a)
    }

    pub fn point_x(&self) -> Elem<R> {
        let mut r = Elem::zero();
        r.inner = self.0[0].clone();
        r
    }

    pub fn point_y(&self) -> Elem<R> {
        let mut r = Elem::zero();
        r.inner = self.0[1].clone();
        r
    }

    pub fn point_z(&self) -> Elem<R> {
        let mut r = Elem::zero();
        r.inner = self.0[2].clone();
        r
    }

    pub fn to_bns(&self) -> [BigUint; 3] {
        self.0.clone()
    }
}

#[derive(Copy, Clone)]
pub struct PublicKey {
    bytes: [u8; PUBLIC_KEY_LEN],
}

impl PublicKey {
    pub fn new(x: &[u8; BN_LENGTH], y: &[u8; BN_LENGTH]) -> Self {
        let mut public = PublicKey {
            bytes: [0; PUBLIC_KEY_LEN],
        };
        public.bytes[0] = 4;
        public.bytes[1..1 + BN_LENGTH].copy_from_slice(x);
        public.bytes[1 + BN_LENGTH..].copy_from_slice(y);

        public
    }

    pub fn bytes_less_safe(&self) -> &[u8] {
        &self.bytes
    }

    pub fn to_point(&self, cctx: &CurveCtx) -> Point {
        let x = BigUint::from_bytes_be(&self.bytes[1..BN_LENGTH + 1]);
        let x_aff = bn_to_mont(&x, cctx);

        let y = BigUint::from_bytes_be(&self.bytes[BN_LENGTH + 1..]);
        let y_aff = bn_to_mont(&y, cctx);

        Point(bn_to_jacobi(&[x_aff, y_aff], cctx))
    }
}

pub const BN_LENGTH: usize = 32;
pub const PUBLIC_KEY_LEN: usize = 1 + (2 * BN_LENGTH);
