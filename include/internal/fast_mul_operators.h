// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "defs.h"

EXTERNC_BEGIN

typedef struct mul_op_s {
  __uint128_t op;
  __uint128_t con;
} mul_op_t;

static inline uint64_t reduce_2q_to_q(const uint64_t val, const uint64_t q)
{
  return (val < q) ? val : val - q;
}

static inline uint64_t reduce_4q_to_2q(const uint64_t val, const uint64_t q)
{
  return (val < 2 * q) ? val : val - 2 * q;
}

static inline uint64_t reduce_4q_to_q(const uint64_t val, const uint64_t q)
{
  return reduce_2q_to_q(reduce_4q_to_2q(val, q), q);
}

static inline uint64_t reduce_8q_to_4q(const uint64_t val, const uint64_t q)
{
  return (val < 4 * q) ? val : val - 4 * q;
}

static inline uint64_t reduce_8q_to_2q(const uint64_t val, const uint64_t q)
{
  return reduce_4q_to_2q(reduce_8q_to_4q(val, q), q);
}

static inline uint64_t reduce_8q_to_q(const uint64_t val, const uint64_t q)
{
  return reduce_2q_to_q(reduce_8q_to_2q(val, q), q);
}

#ifndef L_HIGH_WORD
#  define L_HIGH_WORD HIGH_WORD
#endif

static inline uint64_t
fast_mul_mod_q2(const mul_op_t w, const uint64_t t, const uint64_t q)
{
  const uint64_t Q = L_HIGH_WORD(w.con * t);
  return w.op * t - Q * q;
}

static inline uint64_t
fast_mul_mod_q(const mul_op_t w, const uint64_t t, const uint64_t q)
{
  return reduce_2q_to_q(fast_mul_mod_q2(w, t, q), q);
}

static inline uint64_t fast_dbl_mul_mod_q2(const mul_op_t w1,
                                           const mul_op_t w2,
                                           const uint64_t t1,
                                           const uint64_t t2,
                                           const uint64_t q)
{
  const uint64_t Q = L_HIGH_WORD(w1.con * t1 + w2.con * t2);
  return (t1 * w1.op) + (t2 * w2.op) - (Q * q);
}

static inline void
harvey_fwd_butterfly(uint64_t *X, uint64_t *Y, const mul_op_t w, const uint64_t q)
{
  const uint64_t q2 = q << 1;
  const uint64_t X1 = reduce_4q_to_2q(*X, q);
  const uint64_t T  = fast_mul_mod_q2(w, *Y, q);

  *X = X1 + T;
  *Y = X1 - T + q2;
}

static inline void
harvey_bkw_butterfly(uint64_t *X, uint64_t *Y, const mul_op_t w, const uint64_t q)
{
  const uint64_t q2 = q << 1;
  const uint64_t X1 = reduce_4q_to_2q(*X + *Y, q);
  const uint64_t T  = *X - *Y + q2;

  *X = X1;
  *Y = fast_mul_mod_q2(w, T, q);
}

static inline void harvey_bkw_butterfly_final(uint64_t *     X,
                                              uint64_t *     Y,
                                              const mul_op_t w,
                                              const mul_op_t n_inv,
                                              const uint64_t q)
{
  const uint64_t q2 = q << 1;
  const uint64_t X1 = *X + *Y;
  const uint64_t T  = *X - *Y + q2;

  *X = fast_mul_mod_q(n_inv, X1, q);
  *Y = fast_mul_mod_q(w, T, q);
}

static inline void radix4_fwd_butterfly(uint64_t *     X,
                                        uint64_t *     Y,
                                        uint64_t *     Z,
                                        uint64_t *     T,
                                        const mul_op_t w[5],
                                        const uint64_t q)
{
  const uint64_t q2 = q << 1;
  const uint64_t q4 = q << 2;

  const uint64_t Y1 = fast_dbl_mul_mod_q2(w[1], w[2], *Y, *T, q);
  const uint64_t Y2 = fast_dbl_mul_mod_q2(w[3], w[4], *Y, *T, q);

  const uint64_t T1 = reduce_8q_to_4q(*X, q);
  const uint64_t T2 = fast_mul_mod_q2(w[0], *Z, q);

  *X = (T1 + T2 + Y1);
  *Y = (T1 + T2 - Y1) + q2;
  *Z = (T1 - T2 + Y2) + q2;
  *T = (T1 - T2 - Y2) + q4;
}

static inline void radix4_inv_butterfly(uint64_t *     X,
                                        uint64_t *     Y,
                                        uint64_t *     Z,
                                        uint64_t *     T,
                                        const mul_op_t w[5],
                                        const uint64_t q)
{
  const uint64_t q4 = q << 2;

  const uint64_t T0 = *Z + *T;
  const uint64_t T1 = *X + *Y;

  const uint64_t T2 = q4 + *X - *Y;
  const uint64_t T3 = q4 + *Z - *T;

  *X = reduce_8q_to_2q(T1 + T0, q);
  *Z = fast_mul_mod_q(w[0], q4 + T1 - T0, q);
  *Y = fast_dbl_mul_mod_q2(w[1], w[3], T2, T3, q);
  *T = fast_dbl_mul_mod_q2(w[2], w[4], T2, T3, q);
}

EXTERNC_END
