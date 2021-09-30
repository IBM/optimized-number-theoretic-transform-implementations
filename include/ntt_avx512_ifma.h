// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "defs.h"

EXTERNC_BEGIN

#ifdef AVX512_IFMA_SUPPORT

#  include "avx512.h"

void fwd_ntt_radix4_avx512_ifma_lazy(uint64_t       a[],
                                     uint64_t       N,
                                     uint64_t       q,
                                     const uint64_t w[],
                                     const uint64_t w_con[]);

// Assumption N % 2^6 = 0
static inline void
final_reduce_q8(uint64_t a[], const uint64_t N, const uint64_t q_64)
{
  const __m512i q  = SET1(q_64);
  const __m512i q2 = SET1(q_64 << 1);
  const __m512i q4 = SET1(q_64 << 2);

  // Final reduction
  for(size_t i = 0; i < N; i += 8 * 8) {
    LOOP_UNROLL_8
    for(size_t j = 0; j < 8; j++) {
      __m512i T = reduce_if_greater(LOAD(&a[i + 8 * j]), q4);
      T         = reduce_if_greater(T, q2);
      T         = reduce_if_greater(T, q);
      STORE(&a[i + 8 * j], T);
    }
  }
}

// Assumption N % 2^6 = 0
static inline void
final_reduce_q4(uint64_t a[], const uint64_t N, const uint64_t q_64)
{
  const __m512i q  = SET1(q_64);
  const __m512i q2 = SET1(q_64 << 1);
  // Final reduction
  for(size_t i = 0; i < N; i += 8 * 8) {
    LOOP_UNROLL_8
    for(size_t j = 0; j < 8; j++) {
      __m512i T = reduce_if_greater(LOAD(&a[i + 8 * j]), q2);
      STORE(&a[i + 8 * j], reduce_if_greater(T, q));
    }
  }
}

static inline void fwd_ntt_radix4_avx512_ifma(uint64_t       a[],
                                              const uint64_t N,
                                              const uint64_t q,
                                              const uint64_t w[],
                                              const uint64_t w_con[])
{
  fwd_ntt_radix4_avx512_ifma_lazy(a, N, q, w, w_con);
  final_reduce_q8(a, N, q);
}

void fwd_ntt_r4r2_avx512_ifma_lazy(uint64_t       a[],
                                   uint64_t       N,
                                   uint64_t       q,
                                   const uint64_t w[],
                                   const uint64_t w_con[]);

static inline void fwd_ntt_r4r2_avx512_ifma(uint64_t       a[],
                                            const uint64_t N,
                                            const uint64_t q,
                                            const uint64_t w[],
                                            const uint64_t w_con[])
{
  fwd_ntt_r4r2_avx512_ifma_lazy(a, N, q, w, w_con);
  final_reduce_q4(a, N, q);
}

void fwd_ntt_radix4_avx512_ifma_lazy_unordered(uint64_t       a[],
                                               uint64_t       N,
                                               uint64_t       q,
                                               const uint64_t w[],
                                               const uint64_t w_con[]);

static inline void fwd_ntt_radix4_avx512_ifma_unordered(uint64_t       a[],
                                                        const uint64_t N,
                                                        const uint64_t q,
                                                        const uint64_t w[],
                                                        const uint64_t w_con[])
{
  fwd_ntt_radix4_avx512_ifma_lazy_unordered(a, N, q, w, w_con);
  final_reduce_q8(a, N, q);
}

void fwd_ntt_r2_16_avx512_ifma_lazy(uint64_t       a[],
                                    uint64_t       N,
                                    uint64_t       q,
                                    const uint64_t w[],
                                    const uint64_t w_con[]);

static inline void fwd_ntt_r2_16_avx512_ifma(uint64_t       a[],
                                             const uint64_t N,
                                             const uint64_t q,
                                             const uint64_t w[],
                                             const uint64_t w_con[])
{
  fwd_ntt_r2_16_avx512_ifma_lazy(a, N, q, w, w_con);
  final_reduce_q4(a, N, q);
}

#endif

EXTERNC_END
