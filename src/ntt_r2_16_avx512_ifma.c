// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifdef AVX512_IFMA_SUPPORT

#  include "ntt_avx512_ifma.h"

static inline void fwd16_r2(uint64_t *      a,
                            const uint64_t  m,
                            const uint64_t *w,
                            const uint64_t *w_con,
                            const uint64_t  q_64)
{
  for(size_t j = m; j > 0; --j) {
    size_t              i      = j - 1;
    const uint64_t      w_idx2 = (8 * i) + m;
    const uint64_t      w_idx3 = w_idx2 + (8 * m);
    const uint64_t      w_idx4 = w_idx2 + (16 * m);
    const mul_op_m512_t w1     = {SET1(w[i]), SET1(w_con[i])};
    const mul_op_m512_t w2     = {LOAD(&w[w_idx2]), LOAD(&w_con[w_idx2])};
    const mul_op_m512_t w3     = {LOAD(&w[w_idx3]), LOAD(&w_con[w_idx3])};
    const mul_op_m512_t w4     = {LOAD(&w[w_idx4]), LOAD(&w_con[w_idx4])};

    __m512i X = LOAD(&a[16 * i]);
    __m512i Y = LOAD(&a[16 * i + 8]);
    __m512i T;

    fwd_radix2_butterfly_m512(&X, &Y, &w1, q_64);

    T = SHUF(X, Y, 0x44); // (0, 1 ,2, 3, 80, 81, 82, 83)
    Y = SHUF(X, Y, 0xee); // (4, 5 ,6, 7, 84, 85, 86, 87)
    X = T;

    fwd_radix2_butterfly_m512(&X, &Y, &w2, q_64);

    T = SHUF(X, Y, 0x88); // (0, 1 ,80, 81, 4, 5, 84, 85)
    Y = SHUF(X, Y, 0xdd); // (2, 3, 82, 83, 6, 7, 86, 87)
    X = T;

    fwd_radix2_butterfly_m512(&X, &Y, &w3, q_64);

    __m512i idx1 = SETR(0, 2, 8 + 0, 8 + 2, 4, 6, 8 + 4, 8 + 6);
    __m512i idx2 = SETR(1, 3, 8 + 1, 8 + 3, 5, 7, 8 + 5, 8 + 7);

    T = PERM(X, idx1, Y); // (0, 80 ,2, 82, 4, 84, 6, 86)
    Y = PERM(X, idx2, Y); // (1, 81 ,3, 83, 5, 85, 7, 87)
    X = T;

    fwd_radix2_butterfly_m512(&X, &Y, &w4, q_64);

    STORE(&a[16 * i], UNPACKLO(X, Y));
    STORE(&a[16 * i + 8], UNPACKHI(X, Y));
  }
}

static inline void fwd8_r2(uint64_t *           X_64,
                           uint64_t *           Y_64,
                           const mul_op_m512_t *w,
                           const uint64_t       q_64)
{
  __m512i X = LOAD(X_64);
  __m512i Y = LOAD(Y_64);

  fwd_radix2_butterfly_m512(&X, &Y, w, q_64);

  STORE(X_64, X);
  STORE(Y_64, Y);
}

void fwd_ntt_r2_16_avx512_ifma_lazy(uint64_t       a[],
                                    uint64_t       N,
                                    uint64_t       q,
                                    const uint64_t w[],
                                    const uint64_t w_con[])
{
  size_t m = 1;
  size_t t = N >> 1;

  for(; m < (N >> 4); m <<= 1, t >>= 1) {
    for(size_t j = 0; j < m; j++) {

      const uint64_t      k  = 2 * t * j;
      const mul_op_m512_t w1 = {SET1(w[m + j]), SET1(w_con[m + j])};

      for(size_t i = k; i < k + t; i += 8) {
        fwd8_r2(&a[i], &a[i + t], &w1, q);
      }
    }
  }

  fwd16_r2(a, m, &w[m], &w_con[m], q);
}

#endif
