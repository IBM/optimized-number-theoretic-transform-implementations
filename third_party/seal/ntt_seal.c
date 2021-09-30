// Copyright IBM.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The code belwo was taken and modified from
// https://github.com/microsoft/SEAL/blob/d045f1beff96dff0fccc7fa0c5acb1493a65338c/native/src/seal/util/dwthandler.h
// The license for most of this file is therefore 
//
// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "ntt_seal.h"
#include "defs.h"

// SEAL uses the following Arithmetic inline functions for the NTT implementation

static inline uint64_t add(const uint64_t a, const uint64_t b) { return a + b; }

static inline uint64_t
sub(const uint64_t a, const uint64_t b, const uint64_t two_times_modulus_)
{
  return a + two_times_modulus_ - b;
}

static inline uint64_t mul_root(const uint64_t a,
                                const uint64_t q,
                                const uint64_t w,
                                const uint64_t w_con)
{
  unsigned long long tmp1;
  tmp1 = (unsigned long long)(((__uint128_t)a * (__uint128_t)w_con) >> WORD_SIZE);
  return w * a - tmp1 * q;
}

static inline uint64_t mul_scalar(const uint64_t a,
                                  const uint64_t q,
                                  const uint64_t s,
                                  const uint64_t s_con)
{
  return mul_root(a, q, s, s_con);
}

static inline uint64_t guard(const uint64_t a, const uint64_t two_times_modulus_)
{
  return (a >= two_times_modulus_ ? a - two_times_modulus_ : a);
}

void fwd_ntt_seal_lazy(uint64_t       a[],
                       const uint64_t N,
                       const uint64_t q,
                       const uint64_t w[],
                       const uint64_t w_con[])
{
  // constant transform size
  // Original line: size_t n = size_t(1) << log_n;
  size_t n = N;
  // registers to hold temporary values
  uint64_t u;
  uint64_t v;
  // pointers for faster indexing
  uint64_t *x = NULL;
  uint64_t *y = NULL;
  // variables for indexing
  size_t   gap                = n >> 1;
  size_t   m                  = 1;
  uint64_t two_times_modulus_ = q << 1;

  for(; m < (n >> 1); m <<= 1) {
    size_t offset = 0;
    if(gap < 4) {
      for(size_t i = 0; i < m; i++) {
        ++w;
        ++w_con;
        x = a + offset;
        y = x + gap;
        for(size_t j = 0; j < gap; j++) {
          u    = guard(*x, two_times_modulus_);
          v    = mul_root(*y, q, *w, *w_con);
          *x++ = add(u, v);
          *y++ = sub(u, v, two_times_modulus_);
        }
        offset += gap << 1;
      }
    } else {
      for(size_t i = 0; i < m; i++) {
        ++w;
        ++w_con;
        x = a + offset;
        y = x + gap;
        for(size_t j = 0; j < gap; j += 4) {
          u    = guard(*x, two_times_modulus_);
          v    = mul_root(*y, q, *w, *w_con);
          *x++ = add(u, v);
          *y++ = sub(u, v, two_times_modulus_);

          u    = guard(*x, two_times_modulus_);
          v    = mul_root(*y, q, *w, *w_con);
          *x++ = add(u, v);
          *y++ = sub(u, v, two_times_modulus_);

          u    = guard(*x, two_times_modulus_);
          v    = mul_root(*y, q, *w, *w_con);
          *x++ = add(u, v);
          *y++ = sub(u, v, two_times_modulus_);

          u    = guard(*x, two_times_modulus_);
          v    = mul_root(*y, q, *w, *w_con);
          *x++ = add(u, v);
          *y++ = sub(u, v, two_times_modulus_);
        }
        offset += gap << 1;
      }
    }
    gap >>= 1;
  }

  for(size_t i = 0; i < m; i++) {
    ++w;
    ++w_con;
    u    = guard(a[0], two_times_modulus_);
    v    = mul_root(a[1], q, *w, *w_con);
    a[0] = add(u, v);
    a[1] = sub(u, v, two_times_modulus_);
    a += 2;
  }
}

void inv_ntt_seal(uint64_t       a[],
                  const uint64_t N,
                  const uint64_t q,
                  const uint64_t n_inv,
                  const uint64_t n_inv_con,
                  const uint64_t w[],
                  const uint64_t w_con[])
{
  // constant transform size
  // Original line: size_t n = size_t(1) << log_n;
  size_t n = N;
  // registers to hold temporary values
  uint64_t u;
  uint64_t v;
  // pointers for faster indexing
  uint64_t *x = NULL;
  uint64_t *y = NULL;
  // variables for indexing
  size_t   gap                = 1;
  size_t   m                  = n >> 1;
  uint64_t two_times_modulus_ = q << 1;

  for(; m > 1; m >>= 1) {
    size_t offset = 0;
    if(gap < 4) {
      for(size_t i = 0; i < m; i++) {
        x = a + offset;
        y = x + gap;
        for(size_t j = 0; j < gap; j++) {
          u    = *x;
          v    = *y;
          *x++ = guard(add(u, v), two_times_modulus_);
          *y++ =
            mul_root(sub(u, v, two_times_modulus_), q, w[m + i], w_con[m + i]);
        }
        offset += gap << 1;
      }
    } else {
      for(size_t i = 0; i < m; i++) {
        x = a + offset;
        y = x + gap;
        for(size_t j = 0; j < gap; j += 4) {
          u    = *x;
          v    = *y;
          *x++ = guard(add(u, v), two_times_modulus_);
          *y++ =
            mul_root(sub(u, v, two_times_modulus_), q, w[m + i], w_con[m + i]);

          u    = *x;
          v    = *y;
          *x++ = guard(add(u, v), two_times_modulus_);
          *y++ =
            mul_root(sub(u, v, two_times_modulus_), q, w[m + i], w_con[m + i]);

          u    = *x;
          v    = *y;
          *x++ = guard(add(u, v), two_times_modulus_);
          *y++ =
            mul_root(sub(u, v, two_times_modulus_), q, w[m + i], w_con[m + i]);

          u    = *x;
          v    = *y;
          *x++ = guard(add(u, v), two_times_modulus_);
          *y++ =
            mul_root(sub(u, v, two_times_modulus_), q, w[m + i], w_con[m + i]);
        }
        offset += gap << 1;
      }
    }
    gap <<= 1;
  }

  // Adaption to meet the current code package style
  uint64_t scaled_r     = mul_root(w[1], q, n_inv, n_inv_con);
  uint64_t scaled_r_con = ((__uint128_t)scaled_r << WORD_SIZE) / q;

  x = a;
  y = x + gap;
  if(gap < 4) {
    for(size_t j = 0; j < gap; j++) {
      u = guard(*x, two_times_modulus_);
      v = *y;
      *x++ =
        mul_scalar(guard(add(u, v), two_times_modulus_), q, n_inv, n_inv_con);
      *y++ = mul_root(sub(u, v, two_times_modulus_), q, scaled_r, scaled_r_con);
    }
  } else {
    for(size_t j = 0; j < gap; j += 4) {
      u = guard(*x, two_times_modulus_);
      v = *y;
      *x++ =
        mul_scalar(guard(add(u, v), two_times_modulus_), q, n_inv, n_inv_con);
      *y++ = mul_root(sub(u, v, two_times_modulus_), q, scaled_r, scaled_r_con);

      u = guard(*x, two_times_modulus_);
      v = *y;
      *x++ =
        mul_scalar(guard(add(u, v), two_times_modulus_), q, n_inv, n_inv_con);
      *y++ = mul_root(sub(u, v, two_times_modulus_), q, scaled_r, scaled_r_con);

      u = guard(*x, two_times_modulus_);
      v = *y;
      *x++ =
        mul_scalar(guard(add(u, v), two_times_modulus_), q, n_inv, n_inv_con);
      *y++ = mul_root(sub(u, v, two_times_modulus_), q, scaled_r, scaled_r_con);

      u = guard(*x, two_times_modulus_);
      v = *y;
      *x++ =
        mul_scalar(guard(add(u, v), two_times_modulus_), q, n_inv, n_inv_con);
      *y++ = mul_root(sub(u, v, two_times_modulus_), q, scaled_r, scaled_r_con);
    }
  }

  for(size_t i = 0; i < N; i++) {
    a[i] = (a[i] < q) ? a[i] : a[i] - q;
  }
}
