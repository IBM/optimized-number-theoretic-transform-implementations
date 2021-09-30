// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifdef AVX512_IFMA_SUPPORT

#  include "ntt_avx512_ifma.h"
#  include "ntt_hexl.h"

static inline void _fwd8_r2(uint64_t *           a,
                            __m512i *            X,
                            __m512i *            Y,
                            const mul_op_m512_t *w2,
                            const mul_op_m512_t *w3,
                            const mul_op_m512_t *w4,
                            const uint64_t       q_64)
{
  __m512i T;
  T  = SHUF(*X, *Y, 0x44); // (0, 1 ,2, 3, 80, 81, 82, 83)
  *Y = SHUF(*X, *Y, 0xee); // (4, 5 ,6, 7, 84, 85, 86, 87)
  *X = T;

  fwd_radix2_butterfly_m512(X, Y, w2, q_64);

  T  = SHUF(*X, *Y, 0x88); // (0, 1 ,80, 81, 4, 5, 84, 85)
  *Y = SHUF(*X, *Y, 0xdd); // (2, 3, 82, 83, 6, 7, 86, 87)
  *X = T;

  fwd_radix2_butterfly_m512(X, Y, w3, q_64);

  __m512i idx1 = SETR(0, 2, 8 + 0, 8 + 2, 4, 6, 8 + 4, 8 + 6);
  __m512i idx2 = SETR(1, 3, 8 + 1, 8 + 3, 5, 7, 8 + 5, 8 + 7);

  T  = PERM(*X, idx1, *Y); // (0, 80 ,2, 82, 4, 84, 6, 86)
  *Y = PERM(*X, idx2, *Y); // (1, 81 ,3, 83, 5, 85, 7, 87)
  *X = T;

  fwd_radix2_butterfly_m512(X, Y, w4, q_64);

  STORE(&a[0], UNPACKLO(*X, *Y));
  STORE(&a[8], UNPACKHI(*X, *Y));
}

static inline void fwd8_r2(uint64_t *      a,
                           const uint64_t  m,
                           const uint64_t *w,
                           const uint64_t *w_con,
                           const uint64_t  q_64)
{
  const __m512i q4 = SET1(q_64 << 2);

  LOOP_UNROLL_4
  for(size_t i = 0; i < m; ++i) {
    const uint64_t      w_idx2 = (8 * i);
    const uint64_t      w_idx3 = w_idx2 + (8 * m);
    const uint64_t      w_idx4 = w_idx2 + (16 * m);
    const mul_op_m512_t w2     = {LOAD(&w[w_idx2]), LOAD(&w_con[w_idx2])};
    const mul_op_m512_t w3     = {LOAD(&w[w_idx3]), LOAD(&w_con[w_idx3])};
    const mul_op_m512_t w4     = {LOAD(&w[w_idx4]), LOAD(&w_con[w_idx4])};

    // Radix-4 butterfly leaves values in [0, 8q)
    __m512i X =
      reduce_if_greater(LOAD(&a[16 * i]) & AVX512_IFMA_WORD_SIZE_MASK, q4);
    __m512i Y =
      reduce_if_greater(LOAD(&a[16 * i + 8]) & AVX512_IFMA_WORD_SIZE_MASK, q4);

    _fwd8_r2(&a[16 * i], &X, &Y, &w2, &w3, &w4, q_64);
  }
}

static inline void fwd16_r2(uint64_t *      a,
                            const uint64_t  m,
                            const uint64_t *w,
                            const uint64_t *w_con,
                            const uint64_t  q_64)
{
  const __m512i q4 = SET1(q_64 << 2);

  LOOP_UNROLL_4
  for(size_t j = m; j > 0; --j) {
    size_t              i      = j - 1;
    const uint64_t      w_idx2 = (8 * i) + m;
    const uint64_t      w_idx3 = w_idx2 + (8 * m);
    const uint64_t      w_idx4 = w_idx2 + (16 * m);
    const mul_op_m512_t w1     = {SET1(w[i]), SET1(w_con[i])};
    const mul_op_m512_t w2     = {LOAD(&w[w_idx2]), LOAD(&w_con[w_idx2])};
    const mul_op_m512_t w3     = {LOAD(&w[w_idx3]), LOAD(&w_con[w_idx3])};
    const mul_op_m512_t w4     = {LOAD(&w[w_idx4]), LOAD(&w_con[w_idx4])};

    // Radix-4 butterfly leaves values in [0, 8q)
    __m512i X =
      reduce_if_greater(LOAD(&a[16 * i]) & AVX512_IFMA_WORD_SIZE_MASK, q4);
    __m512i Y =
      reduce_if_greater(LOAD(&a[16 * i + 8]) & AVX512_IFMA_WORD_SIZE_MASK, q4);

    fwd_radix2_butterfly_m512(&X, &Y, &w1, q_64);

    _fwd8_r2(&a[16 * i], &X, &Y, &w2, &w3, &w4, q_64);
  }
}

static inline void collect_roots_fwd8_r4(mul_op_m512_t  w1[5],
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

static inline void fwd8_r4(uint64_t *          X_64,
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

void fwd_ntt_r4r2_avx512_ifma_lazy(uint64_t       a[],
                                   uint64_t       N,
                                   uint64_t       q,
                                   const uint64_t w[],
                                   const uint64_t w_con[])
{
  mul_op_m512_t roots[5];
  size_t        t   = N >> 2;
  size_t        m   = 1;
  size_t        idx = 1;

  for(; t > 4; m <<= 2) {
    for(size_t j = 0; j < m; j++) {
      const uint64_t k = 4 * t * j;
      collect_roots_fwd8_r4(roots, w, w_con, &idx);
      for(size_t i = k; i < k + t; i += 8) {
        fwd8_r4(&a[i], &a[i + t], &a[i + 2 * t], &a[i + 3 * t], roots, q);
      }
    }
    t >>= 2;
  }

  if(HAS_AN_EVEN_POWER(N)) {
    fwd16_r2(a, m, &w[idx], &w_con[idx], q);
  } else {
    m >>= 1;
    fwd8_r2(a, m, &w[idx], &w_con[idx], q);
  }
}

#endif
