// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifdef AVX512_IFMA_SUPPORT

#  include "ntt_avx512_ifma.h"

static inline void collect_roots_fwd1(mul_op_m512_t  w1[5],
                                      const uint64_t w[],
                                      const uint64_t w_con[],
                                      size_t *       idx)
{
  w1[0].op = LOADA(&w[*idx]);
  w1[1].op = LOADA(&w[*idx + 8]);
  w1[2].op = LOADA(&w[*idx + 16]);
  w1[3].op = LOADA(&w[*idx + 24]);
  w1[4].op = LOADA(&w[*idx + 32]);

  w1[0].con = LOADA(&w_con[*idx]);
  w1[1].con = LOADA(&w_con[*idx + 8]);
  w1[2].con = LOADA(&w_con[*idx + 16]);
  w1[3].con = LOADA(&w_con[*idx + 24]);
  w1[4].con = LOADA(&w_con[*idx + 32]);

  *idx += 5 * 8;
}

static inline void collect_roots_fwd4(mul_op_m512_t  w1[5],
                                      const uint64_t w[],
                                      const uint64_t w_con[],
                                      size_t *       idx)
{
  w1[0].op = BROADCAST2HALVES(w[*idx + 0], w[*idx + 1]);
  w1[1].op = BROADCAST2HALVES(w[*idx + 2], w[*idx + 3]);
  w1[2].op = BROADCAST2HALVES(w[*idx + 4], w[*idx + 5]);
  w1[3].op = BROADCAST2HALVES(w[*idx + 6], w[*idx + 7]);
  w1[4].op = BROADCAST2HALVES(w[*idx + 8], w[*idx + 9]);

  w1[0].con = BROADCAST2HALVES(w_con[*idx + 0], w_con[*idx + 1]);
  w1[1].con = BROADCAST2HALVES(w_con[*idx + 2], w_con[*idx + 3]);
  w1[2].con = BROADCAST2HALVES(w_con[*idx + 4], w_con[*idx + 5]);
  w1[3].con = BROADCAST2HALVES(w_con[*idx + 6], w_con[*idx + 7]);
  w1[4].con = BROADCAST2HALVES(w_con[*idx + 8], w_con[*idx + 9]);

  *idx += 10;
}

static inline void collect_roots_fwd8(mul_op_m512_t  w1[5],
                                      const uint64_t w[],
                                      const uint64_t w_con[],
                                      size_t *       idx)
{
  w1[0].op = SET1(w[*idx]);
  w1[1].op = SET1(w[*idx + 1]);
  w1[2].op = SET1(w[*idx + 2]);
  w1[3].op = SET1(w[*idx + 3]);
  w1[4].op = SET1(w[*idx + 4]);

  w1[0].con = SET1(w_con[*idx]);
  w1[1].con = SET1(w_con[*idx + 1]);
  w1[2].con = SET1(w_con[*idx + 2]);
  w1[3].con = SET1(w_con[*idx + 3]);
  w1[4].con = SET1(w_con[*idx + 4]);

  *idx += 5;
}

static inline void
fwd1(uint64_t *a, const mul_op_m512_t w[5], const uint64_t q_64)
{
  const __m512i idx = _mm512_setr_epi64(0, 4, 8, 12, 16, 20, 24, 28);

  __m512i X = _mm512_i64gather_epi64(idx, &a[0 + 0], 8);
  __m512i Y = _mm512_i64gather_epi64(idx, &a[0 + 1], 8);
  __m512i Z = _mm512_i64gather_epi64(idx, &a[0 + 2], 8);
  __m512i T = _mm512_i64gather_epi64(idx, &a[0 + 3], 8);

  fwd_radix4_butterfly_m512(&X, &Y, &Z, &T, w, q_64);

  STORE(&a[0], X);
  STORE(&a[8], Y);
  STORE(&a[16], Z);
  STORE(&a[24], T);
}

static inline void
fwd4(uint64_t *a, const mul_op_m512_t w[5], const uint64_t q_64)
{
  __m512i X1 = LOAD(&a[0]);
  __m512i Y1 = LOAD(&a[8]);
  __m512i Z1 = LOAD(&a[16]);
  __m512i T1 = LOAD(&a[24]);

  __m512i X = _mm512_shuffle_i64x2(X1, Z1, 0x44);
  __m512i Y = _mm512_shuffle_i64x2(X1, Z1, 0xee);
  __m512i Z = _mm512_shuffle_i64x2(Y1, T1, 0x44);
  __m512i T = _mm512_shuffle_i64x2(Y1, T1, 0xee);

  fwd_radix4_butterfly_m512(&X, &Y, &Z, &T, w, q_64);

  STORE(&a[0], X);
  STORE(&a[8], Y);
  STORE(&a[16], Z);
  STORE(&a[24], T);
}

static inline void fwd8(uint64_t *          X_64,
                        uint64_t *          Y_64,
                        uint64_t *          Z_64,
                        uint64_t *          T_64,
                        const mul_op_m512_t w[5],
                        const uint64_t      q_64)
{
  __m512i X = LOAD(X_64);
  __m512i Y = LOAD(Y_64);
  __m512i Z = LOAD(Z_64);
  __m512i T = LOAD(T_64);

  fwd_radix4_butterfly_m512(&X, &Y, &Z, &T, w, q_64);

  STORE(X_64, X);
  STORE(Y_64, Y);
  STORE(Z_64, Z);
  STORE(T_64, T);
}

void fwd_ntt_radix4_avx512_ifma_lazy_unordered(uint64_t       a[],
                                               const uint64_t N,
                                               const uint64_t q,
                                               const uint64_t w[],
                                               const uint64_t w_con[])
{
  mul_op_m512_t roots[5];
  size_t        bound_r4 = N;
  size_t        t        = N >> 1;
  size_t        m        = 1;
  size_t        idx      = 1;

  // Check whether N=2^m where m is odd.
  // If not perform extra radix-2 iteration.
  if(!HAS_AN_EVEN_POWER(N)) {
    const mul_op_m512_t w1 = {SET1(w[1]), SET1(w_con[1])};

    for(size_t j = 0; j < t; j += 8) {
      __m512i X = LOAD(&a[j]);
      __m512i Y = LOAD(&a[j + t]);

      fwd_radix2_butterfly_m512(&X, &Y, &w1, q);

      STORE(&a[j], X);
      STORE(&a[j + t], Y);
    }
    bound_r4 >>= 1;
    t >>= 1;
    m <<= 1;
    idx++;
  }

  // Adjust to radix-4
  t >>= 1;

  for(; m < bound_r4; m <<= 2) {
    if(t >= 8) {
      for(size_t j = 0; j < m; j++) {
        const uint64_t k = 4 * t * j;
        collect_roots_fwd8(roots, w, w_con, &idx);
        for(size_t i = k; i < k + t; i += 8) {
          fwd8(&a[i], &a[i + t], &a[i + 2 * t], &a[i + 3 * t], roots, q);
        }
      }
    } else if(t == 4) {
      for(size_t j = 0; j < m; j += 2) {
        collect_roots_fwd4(roots, w, w_con, &idx);
        fwd4(&a[4 * 4 * j], roots, q);
      }
    } else {
      // Align on an 8-qw boundary
      idx = ((idx >> 3) << 3) + 8;

      for(size_t j = 0; j < m;) {
        LOOP_UNROLL_4
        for(size_t k = 0; k < 4; k += 1, j += 8) {
          collect_roots_fwd1(roots, w, w_con, &idx);
          fwd1(&a[4 * j], roots, q);
        }
      }
    }
    t >>= 2;
  }
}

#endif
