// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "defs.h"

EXTERNC_BEGIN

#include <immintrin.h>

#define ADD(a, b) _mm512_add_epi64(a, b)
#define SUB(a, b) _mm512_sub_epi64(a, b)
#define MIN(a, b) _mm512_min_epu64(a, b)

#define SHUF(a, b, mask) _mm512_shuffle_i64x2((a), (b), (mask))
#define PERM(a, idx, b)  _mm512_permutex2var_epi64((a), (idx), (b))
#define UNPACKLO(a, b)   _mm512_unpacklo_epi64((a), (b))
#define UNPACKHI(a, b)   _mm512_unpackhi_epi64((a), (b))

#define MADDLO(accum, op1, op2) _mm512_madd52lo_epu64((accum), (op1), (op2))
#define MADDHI(accum, op1, op2) _mm512_madd52hi_epu64((accum), (op1), (op2))

#define SET1(val) _mm512_set1_epi64(val)
#define SETR(a, b, c, d, e, f, g, h) \
  _mm512_setr_epi64((a), (b), (c), (d), (e), (f), (g), (h))
#define LOAD(mem)       _mm512_loadu_epi64((mem))
#define STORE(mem, reg) _mm512_storeu_epi64((mem), (reg))

#define GATHER(idx, mem, scale) _mm512_i64gather_epi64((idx), (mem), (scale))
#define SCATTER(mem, idx, reg, scale) \
  _mm512_i64scatter_epi64((mem), (idx), (reg), scale)

#define SET1_256(val)            _mm256_set1_epi64x(val)
#define BROADCAST2HALVES(v1, v2) _mm512_inserti64x4(SET1((v1)), SET1_256((v2)), 1)

typedef struct mul_op_m512_s {
  __m512i op;
  __m512i con;
} mul_op_m512_t;

static inline __m512i reduce_if_greater(const __m512i val, const __m512i mod)
{
  return MIN(val, SUB(val, mod));
}

static inline __m512i
fast_mul_mod_q2_m512(const mul_op_m512_t w, const __m512i t, const __m512i neg_q)
{
  const __m512i zero = SET1(0);
  const __m512i Q    = MADDHI(zero, w.con, t);
  const __m512i tmp  = MADDLO(zero, Q, neg_q);
  return MADDLO(tmp, w.op, t) & AVX512_IFMA_WORD_SIZE_MASK;
}

static inline __m512i fast_dbl_mul_mod_q2_m512(const mul_op_m512_t w1,
                                               const mul_op_m512_t w2,
                                               const __m512i       t1,
                                               const __m512i       t2,
                                               const __m512i       neg_q)
{
  const __m512i zero = SET1(0);
  const __m512i TMP1 = MADDHI(zero, w1.con, t1);
  const __m512i Q    = MADDHI(TMP1, w2.con, t2);

  const __m512i TMP2 = MADDLO(zero, w1.op, t1);
  const __m512i TMP3 = MADDLO(TMP2, w2.op, t2);
  const __m512i TMP4 = MADDLO(TMP3, Q, neg_q);
  return TMP4 & AVX512_IFMA_WORD_SIZE_MASK;
}

static inline void fwd_radix4_butterfly_m512(__m512i *           X,
                                             __m512i *           Y,
                                             __m512i *           Z,
                                             __m512i *           T,
                                             const mul_op_m512_t w[5],
                                             const uint64_t      q_64)
{
  const __m512i neg_q = SET1(-1 * q_64);
  const __m512i q2    = SET1(q_64 << 1);
  const __m512i q4    = SET1(q_64 << 2);

  const __m512i T1 = reduce_if_greater(*X, q4);
  const __m512i T2 = fast_mul_mod_q2_m512(w[0], *Z, neg_q);

  const __m512i Y1 = fast_dbl_mul_mod_q2_m512(w[1], w[2], *Y, *T, neg_q);
  const __m512i Y2 = fast_dbl_mul_mod_q2_m512(w[3], w[4], *Y, *T, neg_q);

  const __m512i T3 = ADD(T1, T2);
  const __m512i T4 = SUB(T1, T2);

  *X = ADD(T3, Y1);
  *Y = ADD(SUB(q2, Y1), T3);
  *Z = ADD(ADD(q2, Y2), T4);
  *T = ADD(SUB(q4, Y2), T4);
}

static inline void fwd_radix2_butterfly_m512(__m512i *            X,
                                             __m512i *            Y,
                                             const mul_op_m512_t *w,
                                             const uint64_t       q_64)
{
  const __m512i neg_q = SET1(-1 * q_64);
  const __m512i q2    = SET1(q_64 << 1);

  *X              = reduce_if_greater(*X, q2);
  const __m512i T = fast_mul_mod_q2_m512(*w, *Y, neg_q);

  *Y = ADD(SUB(q2, T), *X);
  *X = ADD(*X, T);
}

EXTERNC_END
