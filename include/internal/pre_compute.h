// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string.h>

#include "defs.h"

EXTERNC_BEGIN

// We dont care about the performance of the functions
// in this file, as all of their output serves as
// precomputaitons that we can cache for the NTT computations.

static inline uint64_t bit_rev_idx(uint64_t idx, uint64_t width)
{
  uint64_t ret = 0;
  while(width > 0) {
    width--;
    ret += ((idx & 1) << width);
    idx >>= 1;
  }

  return ret;
}

static inline void bit_rev(uint64_t       w_powers[],
                           const uint64_t w[],
                           const uint64_t N,
                           const uint64_t width)
{
  for(size_t i = 0; i < N; i++) {
    w_powers[bit_rev_idx(i, width)] = w[i];
  }
}

static inline void calc_w(uint64_t       w_powers_rev[],
                          const uint64_t w,
                          const uint64_t N,
                          const uint64_t q,
                          const uint64_t width)
{
  uint64_t w_powers[N];
  w_powers[0] = 1;
  for(size_t i = 1; i < N; i++) {
    w_powers[i] = (uint64_t)(((__uint128_t)w_powers[i - 1] * w) % q);
  }

  bit_rev(w_powers_rev, w_powers, N, width);
}

static inline void calc_w_inv(uint64_t       w_inv_rev[],
                              const uint64_t w_inv,
                              const uint64_t N,
                              const uint64_t q,
                              const uint64_t width)
{
  uint64_t w_inv_powers[N];
  w_inv_powers[0] = 1;
  for(size_t i = 1; i < N; i++) {
    w_inv_powers[i] = (uint64_t)(((__uint128_t)w_inv_powers[i - 1] * w_inv) % q);
  }

  bit_rev(w_inv_rev, w_inv_powers, N, width);
}

static inline void calc_w_con(uint64_t       w_con[],
                              const uint64_t w[],
                              const uint64_t N,
                              const uint64_t q,
                              const uint64_t word_size)
{
  for(size_t i = 0; i < N; i++) {
    w_con[i] = ((__uint128_t)w[i] << word_size) / q;
  }
}

static uint64_t
calc_ninv_con(const uint64_t Ninv, const uint64_t q, const uint64_t word_size)
{
  return ((__uint128_t)Ninv << word_size) / q;
}

static inline void expand_w(uint64_t       w_expanded[],
                            const uint64_t w[],
                            const uint64_t N,
                            const uint64_t q)
{
  w_expanded[0] = w[0];
  w_expanded[1] = 0;
  w_expanded[2] = w[1];
  w_expanded[3] = 0;
  for(size_t i = 4; i < 2 * N; i += 2) {
    w_expanded[i] = w[i / 2];

    if(i % 4 == 0) {
      const __uint128_t t = w_expanded[i / 2];
      w_expanded[i + 1]   = (t * w[i / 2]) % q;
    } else {
      const __uint128_t t = w_expanded[(i - 2) / 2];
      w_expanded[i + 1]   = q - ((t * w[i / 2]) % q);
    }
  }
}

#ifdef AVX512_IFMA_SUPPORT

static inline void
expand_w_hexl(uint64_t w_expanded[], const uint64_t w[], const uint64_t N)
{
  size_t idx = 0;

  memcpy(&w_expanded[idx], w, N / 8 * sizeof(uint64_t));
  idx += N / 8;

  // Duplicate four times for FwdT4
  for(size_t i = 0; i < (N / 8); i++) {
    w_expanded[idx]     = w[(N / 8) + i];
    w_expanded[idx + 1] = w[(N / 8) + i];
    w_expanded[idx + 2] = w[(N / 8) + i];
    w_expanded[idx + 3] = w[(N / 8) + i];
    idx += 4;
  }

  // Duplicate four times for FwdT2
  for(size_t i = 0; i < (N / 4); i++) {
    w_expanded[idx]     = w[(N / 4) + i];
    w_expanded[idx + 1] = w[(N / 4) + i];
    idx += 2;
  }

  memcpy(&w_expanded[idx], &w[N / 2], N / 2 * sizeof(uint64_t));
  idx += N / 2;

  memset(&w_expanded[idx], 0, (2 * N - idx) * sizeof(uint64_t));
}

static inline void permute_w(uint64_t in_out[8])
{
  uint64_t t[8];
  memcpy(t, in_out, 8 * sizeof(uint64_t));

  in_out[0] = t[0];
  in_out[1] = t[4];
  in_out[2] = t[1];
  in_out[3] = t[5];
  in_out[4] = t[2];
  in_out[5] = t[6];
  in_out[6] = t[3];
  in_out[7] = t[7];
}

static inline void expand_w_r4_avx512_ifma(uint64_t       w_expanded[],
                                           const uint64_t w[],
                                           const uint64_t N,
                                           const uint64_t q,
                                           const uint64_t unordered)
{
  size_t w_idx     = 1;
  size_t new_w_idx = 1;

  w_expanded[0] = 0;

  // FWD8
  if(HAS_AN_EVEN_POWER(N)) {
    for(size_t m = 1; w_idx < (N >> 5); m <<= 2) {
      for(size_t i = 0; i < m; i++, w_idx++) {
        const __uint128_t w1    = w[w_idx];
        const __uint128_t w2    = w[2 * w_idx];
        const __uint128_t w3    = w[2 * w_idx + 1];
        w_expanded[new_w_idx++] = w1;
        w_expanded[new_w_idx++] = w2;
        w_expanded[new_w_idx++] = (w1 * w2) % q;
        w_expanded[new_w_idx++] = w3;
        w_expanded[new_w_idx++] = q - ((w1 * w3) % q);
      }
      w_idx = 4 * m;
    }
  } else {
    // First radix-2 iteration
    w_expanded[new_w_idx++] = w[w_idx++];

    for(size_t m = 2; w_idx < (N >> 5); m <<= 2) {
      for(size_t i = 0; i < m; i++, w_idx++) {
        const __uint128_t w1    = w[w_idx];
        const __uint128_t w2    = w[2 * w_idx];
        const __uint128_t w3    = w[2 * w_idx + 1];
        w_expanded[new_w_idx++] = w1;
        w_expanded[new_w_idx++] = w2;
        w_expanded[new_w_idx++] = (w1 * w2) % q;
        w_expanded[new_w_idx++] = w3;
        w_expanded[new_w_idx++] = q - ((w1 * w3) % q);
      }
      w_idx = 4 * m;
    }
  }

  // FWD4
  for(w_idx = (N >> 4); w_idx < (N >> 3); w_idx += 2) {
    const uint64_t k        = 2 * w_idx;
    w_expanded[new_w_idx++] = w[w_idx];
    w_expanded[new_w_idx++] = w[w_idx + 1];
    w_expanded[new_w_idx++] = w[k];
    w_expanded[new_w_idx++] = w[k + 2];
    w_expanded[new_w_idx++] = (w[w_idx] * w[k]) % q;
    w_expanded[new_w_idx++] = (w[w_idx + 1] * w[k + 2]) % q;
    w_expanded[new_w_idx++] = w[k + 1];
    w_expanded[new_w_idx++] = w[k + 2 + 1];
    w_expanded[new_w_idx++] = q - ((w[w_idx] * w[k + 1]) % q);
    w_expanded[new_w_idx++] = q - ((w[w_idx + 1] * w[k + 3]) % q);
  }

  // FWD1
  for(w_idx = (N >> 2); w_idx < (N >> 1); w_idx += 8) {
    // W1
    for(size_t i = 0; i < 8; i++) {
      w_expanded[new_w_idx++] = w[w_idx + i];
    }
    // W2
    for(size_t i = 0; i < 8; i++) {
      w_expanded[new_w_idx++] = w[2 * (w_idx + i)];
    }
    // W3
    for(size_t i = 0; i < 8; i++) {
      w_expanded[new_w_idx++] = (w[w_idx + i] * w[2 * (w_idx + i)]) % q;
    }
    // W4
    for(size_t i = 0; i < 8; i++) {
      w_expanded[new_w_idx++] = w[2 * (w_idx + i) + 1];
    }
    // W5
    for(size_t i = 0; i < 8; i++) {
      w_expanded[new_w_idx++] = q - ((w[w_idx + i] * w[2 * (w_idx + i) + 1]) % q);
    }

    // Need to permute values
    if(unordered) {
      permute_w(&w_expanded[new_w_idx - 8 * 5]);
      permute_w(&w_expanded[new_w_idx - 8 * 4]);
      permute_w(&w_expanded[new_w_idx - 8 * 3]);
      permute_w(&w_expanded[new_w_idx - 8 * 2]);
      permute_w(&w_expanded[new_w_idx - 8 * 1]);
    }
  }

  memset(&w_expanded[new_w_idx], 0, ((5 * N) - new_w_idx) * sizeof(uint64_t));
}

static inline void expand_w_r4r2_avx512_ifma(uint64_t       w_expanded[],
                                             const uint64_t w[],
                                             const uint64_t N,
                                             const uint64_t q)
{
  size_t w_idx     = 1;
  size_t new_w_idx = 1;
  size_t t         = N >> 4;

  w_expanded[0] = 0;

  // FWD8 in radix4
  for(size_t m = 1; w_idx < t; m <<= 2) {
    for(size_t i = 0; i < m; i++, w_idx++) {
      const __uint128_t w1    = w[w_idx];
      const __uint128_t w2    = w[2 * w_idx];
      const __uint128_t w3    = w[2 * w_idx + 1];
      w_expanded[new_w_idx++] = w1;
      w_expanded[new_w_idx++] = w2;
      w_expanded[new_w_idx++] = (w1 * w2) % q;
      w_expanded[new_w_idx++] = w3;
      w_expanded[new_w_idx++] = q - ((w1 * w3) % q);
    }
    w_idx = 4 * m;
  }

  if(HAS_AN_EVEN_POWER(N)) {
    // FWD8 in radix2
    memcpy(&w_expanded[new_w_idx], &w[w_idx], t * sizeof(uint64_t));
    new_w_idx += t;
  }

  t <<= 1;

  // Duplicate four times for FwdT4
  for(size_t i = 0; i < t; i++) {
    w_expanded[new_w_idx++] = w[t + i];
    w_expanded[new_w_idx++] = w[t + i];
    w_expanded[new_w_idx++] = w[t + i];
    w_expanded[new_w_idx++] = w[t + i];
  }
  t <<= 1;

  // Duplicate four times for FwdT2
  for(size_t i = 0; i < t; i += 4) {
    w_expanded[new_w_idx++] = w[t + i + 0];
    w_expanded[new_w_idx++] = w[t + i + 0];
    w_expanded[new_w_idx++] = w[t + i + 2];
    w_expanded[new_w_idx++] = w[t + i + 2];
    w_expanded[new_w_idx++] = w[t + i + 1];
    w_expanded[new_w_idx++] = w[t + i + 1];
    w_expanded[new_w_idx++] = w[t + i + 3];
    w_expanded[new_w_idx++] = w[t + i + 3];
  }
  t <<= 1;

  for(size_t i = 0; i < t; i += 8) {
    w_expanded[new_w_idx++] = w[t + i + 0];
    w_expanded[new_w_idx++] = w[t + i + 4];
    w_expanded[new_w_idx++] = w[t + i + 1];
    w_expanded[new_w_idx++] = w[t + i + 5];
    w_expanded[new_w_idx++] = w[t + i + 2];
    w_expanded[new_w_idx++] = w[t + i + 6];
    w_expanded[new_w_idx++] = w[t + i + 3];
    w_expanded[new_w_idx++] = w[t + i + 7];
  }

  memset(&w_expanded[new_w_idx], 0, ((5 * N) - new_w_idx) * sizeof(uint64_t));
}

static inline void expand_w_r2_16_avx512_ifma(uint64_t       w_expanded[],
                                              const uint64_t w[],
                                              const uint64_t N)
{
  size_t t         = N >> 3;
  size_t new_w_idx = t;

  memcpy(w_expanded, w, t * sizeof(uint64_t));

  // Duplicate four times for FwdT4
  for(size_t i = 0; i < t; i++) {
    w_expanded[new_w_idx++] = w[t + i];
    w_expanded[new_w_idx++] = w[t + i];
    w_expanded[new_w_idx++] = w[t + i];
    w_expanded[new_w_idx++] = w[t + i];
  }
  t <<= 1;

  // Duplicate four times for FwdT2
  for(size_t i = 0; i < t; i += 4) {
    w_expanded[new_w_idx++] = w[t + i + 0];
    w_expanded[new_w_idx++] = w[t + i + 0];
    w_expanded[new_w_idx++] = w[t + i + 2];
    w_expanded[new_w_idx++] = w[t + i + 2];
    w_expanded[new_w_idx++] = w[t + i + 1];
    w_expanded[new_w_idx++] = w[t + i + 1];
    w_expanded[new_w_idx++] = w[t + i + 3];
    w_expanded[new_w_idx++] = w[t + i + 3];
  }
  t <<= 1;

  for(size_t i = 0; i < t; i += 8) {
    w_expanded[new_w_idx++] = w[t + i + 0];
    w_expanded[new_w_idx++] = w[t + i + 4];
    w_expanded[new_w_idx++] = w[t + i + 1];
    w_expanded[new_w_idx++] = w[t + i + 5];
    w_expanded[new_w_idx++] = w[t + i + 2];
    w_expanded[new_w_idx++] = w[t + i + 6];
    w_expanded[new_w_idx++] = w[t + i + 3];
    w_expanded[new_w_idx++] = w[t + i + 7];
  }
}

#endif

EXTERNC_END
