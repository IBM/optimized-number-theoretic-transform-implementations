// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <vecintrin.h>

#define L_HIGH_WORD HIGH_VMSL_WORD

#include "fast_mul_operators.h"
#include "ntt_radix4_s390x_vef.h"

#define UL_VMSL_Z(a, b, ctx) (ul_vec) vec_msum_u128(a, b, (ctx)->zero, 0)
#define UL_VMSL              (ul_vec) vec_msum_u128
#define VEC_BYTES            16

typedef vector unsigned long long ul_vec;
typedef vector unsigned char      uc_vec;

typedef struct loop_ctx_s {
  mul_op_t w1;
  ul_vec   r1;
  ul_vec   r2;
  ul_vec   r1_con;
  ul_vec   r2_con;

  // In order not to hit a performance slowdown when using VMSL,
  // we mask -q to be smaller than 2^56. This should be fine in cases
  // where we only consider the lower 56 bits of the VMSL output.
  uint64_t neg_q;

  // q2 = 2*q, q4 = 4*q
  uint64_t q2;
  uint64_t q4;
  ul_vec   q2_vec;
  ul_vec   q4_vec;

  uc_vec zero;
} loop_ctx_t;

/******************************
       Single input
******************************/

static inline ul_vec extended_shoup_multiply(const ul_vec      r,
                                             const ul_vec      r_con,
                                             const ul_vec      YT,
                                             const ul_vec      neg_q,
                                             const loop_ctx_t *ctx)
{
  ul_vec t1 = UL_VMSL_Z(r_con, YT, ctx);
  ul_vec t2 = UL_VMSL_Z(r, YT, ctx);
  t1        = vec_sld((ul_vec)ctx->zero, t1, VEC_BYTES - (VMSL_WORD_SIZE / 8));
  return UL_VMSL(neg_q, t1, (uc_vec)t2, 0);
}

static inline void single_fwd_butterfly(uint64_t          a[],
                                        const uint64_t    q,
                                        const loop_ctx_t *ctx,
                                        const size_t      i,
                                        const uint64_t    t)
{
  const uint64_t X = reduce_8q_to_4q(a[i], q);
  const uint64_t Z = fast_mul_mod_q2(ctx->w1, a[i + 2 * t], q);

  // Load Y and T
  ul_vec YT = {a[i + t], a[i + 3 * t]};

  ul_vec t1         = {0, ctx->neg_q};
  ul_vec r4_1       = extended_shoup_multiply(ctx->r1, ctx->r1_con, YT, t1, ctx);
  ul_vec r5_1       = extended_shoup_multiply(ctx->r2, ctx->r2_con, YT, t1, ctx);
  const uint64_t Y1 = LOW_VMSL_WORD(r4_1[1]);
  const uint64_t Y2 = LOW_VMSL_WORD(r5_1[1]);

  // Store the results
  a[i]         = (X + Z + Y1);
  a[i + t]     = ctx->q2 + (X + Z - Y1);
  a[i + 2 * t] = ctx->q2 + (X - Z + Y2);
  a[i + 3 * t] = ctx->q4 + (X - Z - Y2);
}

// Assumption t is even.
static inline void double_fwd_butterfly(uint64_t          a[],
                                        const loop_ctx_t *ctx,
                                        const size_t      i,
                                        const uint64_t    t)
{
  // We use vector pointers to load two coefficients of "a" at once.
  ul_vec *a_vec = (ul_vec *)&a[i];

  // Load X and Z, note that because we are using vector pointers
  // t is divided by 2.
  ul_vec X = a_vec[0];
  ul_vec Z = a_vec[2 * (t >> 1)];

  // Load Y and T of each of the two iterations into YT1, YT2, respectively.
  // Note that here we use t instead of t>>1 as "a" is a uint64_t.
  ul_vec YT1 = {a[i + t], a[i + 3 * t]};
  ul_vec YT2 = {a[i + t + 1], a[i + 3 * t + 1]};

  // X = X mod 4q
  X = vec_sel(X, X - ctx->q4_vec, vec_cmpge(X, ctx->q4_vec));

  // Simple Shoup multiply on two elements in parallel.
  ul_vec t1  = {ctx->w1.op, ctx->neg_q};
  ul_vec Q_1 = {Z[0], HIGH_VMSL_WORD(ctx->w1.con * Z[0])};
  ul_vec Q_2 = {Z[1], HIGH_VMSL_WORD(ctx->w1.con * Z[1])};
  Q_1        = UL_VMSL_Z(t1, Q_1, ctx);
  Q_2        = UL_VMSL_Z(t1, Q_2, ctx);
  Z          = LOW_VMSL_WORD(vec_mergel(Q_1, Q_2));

  // Extended Shoup multiply on two elements in parallel.
  t1[0]       = 0;
  ul_vec r4_1 = extended_shoup_multiply(ctx->r1, ctx->r1_con, YT1, t1, ctx);
  ul_vec r4_2 = extended_shoup_multiply(ctx->r1, ctx->r1_con, YT2, t1, ctx);
  ul_vec r5_1 = extended_shoup_multiply(ctx->r2, ctx->r2_con, YT1, t1, ctx);
  ul_vec r5_2 = extended_shoup_multiply(ctx->r2, ctx->r2_con, YT2, t1, ctx);

  const ul_vec Y1 = LOW_VMSL_WORD(vec_mergel(r4_1, r4_2));
  const ul_vec Y2 = LOW_VMSL_WORD(vec_mergel(r5_1, r5_2));

  // Store the results
  a_vec[0 * (t >> 1)] = (X + Z + Y1);
  a_vec[1 * (t >> 1)] = ctx->q2_vec + (X + Z - Y1);
  a_vec[2 * (t >> 1)] = ctx->q2_vec + (X - Z + Y2);
  a_vec[3 * (t >> 1)] = ctx->q4_vec + (X - Z - Y2);
}

void fwd_ntt_radix4_intrinsic_lazy(uint64_t       a[],
                                   const uint64_t N,
                                   const uint64_t q,
                                   const uint64_t w[],
                                   const uint64_t w_con[])
{
  const size_t bound_r4 = HAS_AN_EVEN_POWER(N) ? N : (N >> 1);

  for(size_t m = 1, t = (N >> 2); m < bound_r4; m <<= 2, t >>= 2) {
    for(size_t j = 0; j < m; j++) {
      const size_t     m1  = 2 * (m + j);
      const size_t     m2  = 2 * m1;
      const size_t     k   = 4 * t * j;
      const loop_ctx_t ctx = {.w1     = {w[m1], w_con[m1]},
                              .r1     = {w[m2 + 0], w[m2 + 1]},
                              .r2     = {w[m2 + 2], w[m2 + 3]},
                              .r1_con = {w_con[m2 + 0], w_con[m2 + 1]},
                              .r2_con = {w_con[m2 + 2], w_con[m2 + 3]},
                              .zero   = {0},
                              .neg_q  = (-1 * q) & VMSL_WORD_SIZE_MASK,
                              .q2     = 2 * q,
                              .q4     = 4 * q,
                              .q2_vec = {2 * q, 2 * q},
                              .q4_vec = {4 * q, 4 * q}};

      if(t == 1) {
        for(size_t i = k; i < k + t; i++) {
          single_fwd_butterfly(a, q, &ctx, i, t);
        }
      } else {
        for(size_t i = k; i < k + t; i += 2) {
          double_fwd_butterfly(a, &ctx, i, t);
        }
      }
    }
  }

  // Check whether N=2^m where m is odd.
  if(HAS_AN_EVEN_POWER(N)) {
    return;
  }

  for(size_t i = 0; i < N; i += 2) {
    const mul_op_t w1 = {w[i + N], w_con[i + N]};

    a[i] = reduce_8q_to_4q(a[i], q);
    harvey_fwd_butterfly(&a[i], &a[i + 1], w1, q);
  }
}

static inline void single_inv_butterfly(uint64_t          a[],
                                        const uint64_t    q,
                                        const loop_ctx_t *ctx,
                                        const size_t      i,
                                        const uint64_t    t)
{
  const uint64_t X = a[i + 0 * t];
  const uint64_t Y = a[i + 1 * t];
  const uint64_t Z = a[i + 2 * t];
  const uint64_t T = a[i + 3 * t];

  const uint64_t T0 = Z + T;
  const uint64_t T1 = X + Y;
  const uint64_t T4 = fast_mul_mod_q2(ctx->w1, ctx->q4 + T1 - T0, q);

  ul_vec T23  = {ctx->q4 + X - Y, ctx->q4 + Z - T};
  ul_vec t1   = {0, ctx->neg_q};
  ul_vec r4_1 = extended_shoup_multiply(ctx->r1, ctx->r1_con, T23, t1, ctx);
  ul_vec r5_1 = extended_shoup_multiply(ctx->r2, ctx->r2_con, T23, t1, ctx);

  const uint64_t Y1 = LOW_VMSL_WORD(r4_1[1]);
  const uint64_t Y2 = LOW_VMSL_WORD(r5_1[1]);

  a[i + 0 * t] = reduce_8q_to_2q(T1 + T0, q);
  a[i + 1 * t] = Y1;
  a[i + 2 * t] = T4;
  a[i + 3 * t] = Y2;
}

// Assumption t is even.
static inline void double_inv_butterfly(uint64_t          a[],
                                        const loop_ctx_t *ctx,
                                        const size_t      i,
                                        const size_t      t)
{
  // We use vector pointers to load two coefficients of "a" at once.
  ul_vec *a_vec = (ul_vec *)&a[i];

  const ul_vec X = a_vec[0 * t];
  const ul_vec Y = a_vec[1 * t];
  const ul_vec Z = a_vec[2 * t];
  const ul_vec T = a_vec[3 * t];

  const ul_vec T0 = Z + T;
  ul_vec       T1 = X + Y;
  ul_vec       T4 = ctx->q4_vec + T1 - T0;

  // Simple Shoup multiply on two elements in parallel.
  ul_vec t1  = {ctx->w1.op, ctx->neg_q};
  ul_vec Q_1 = {T4[0], HIGH_VMSL_WORD(ctx->w1.con * T4[0])};
  ul_vec Q_2 = {T4[1], HIGH_VMSL_WORD(ctx->w1.con * T4[1])};
  Q_1        = UL_VMSL_Z(t1, Q_1, ctx);
  Q_2        = UL_VMSL_Z(t1, Q_2, ctx);
  T4         = LOW_VMSL_WORD(vec_mergel(Q_1, Q_2));

  ul_vec T2   = ctx->q4_vec + X - Y;
  ul_vec T3   = ctx->q4_vec + Z - T;
  ul_vec T23a = {T2[0], T3[0]};
  ul_vec T23b = {T2[1], T3[1]};

  // Extended Shoup multiply on two elements in parallel.
  t1[0]       = 0;
  ul_vec r4_1 = extended_shoup_multiply(ctx->r1, ctx->r1_con, T23a, t1, ctx);
  ul_vec r4_2 = extended_shoup_multiply(ctx->r1, ctx->r1_con, T23b, t1, ctx);
  ul_vec r5_1 = extended_shoup_multiply(ctx->r2, ctx->r2_con, T23a, t1, ctx);
  ul_vec r5_2 = extended_shoup_multiply(ctx->r2, ctx->r2_con, T23b, t1, ctx);

  const ul_vec Y1 = LOW_VMSL_WORD(vec_mergel(r4_1, r4_2));
  const ul_vec Y2 = LOW_VMSL_WORD(vec_mergel(r5_1, r5_2));

  T1           = T1 + T0;
  T1           = vec_sel(T1, T1 - ctx->q4_vec, vec_cmpge(T1, ctx->q4_vec));
  T1           = vec_sel(T1, T1 - ctx->q2_vec, vec_cmpge(T1, ctx->q2_vec));
  a_vec[0 * t] = T1;
  a_vec[1 * t] = Y1;
  a_vec[2 * t] = T4;
  a_vec[3 * t] = Y2;
}

void inv_ntt_radix4_intrinsic(uint64_t       a[],
                              const uint64_t N,
                              const uint64_t q,
                              const mul_op_t n_inv,
                              const uint64_t w[],
                              const uint64_t w_con[])
{
  size_t t = 1;
  size_t m = N;

  // Check whether N=2^m where m is odd.
  if(HAS_AN_EVEN_POWER(N)) {
    for(size_t i = 0; i < N; i++) {
      a[i] = reduce_8q_to_2q(a[i], q);
    }
  } else {
    // Perform the first iteration as a radix-2 iteration.
    for(size_t i = 0; i < N; i += 2) {
      const mul_op_t w1 = {w[i + N], w_con[i + N]};
      a[i]              = reduce_8q_to_4q(a[i], q);
      harvey_bkw_butterfly(&a[i], &a[i + 1], w1, q);
    }
    m >>= 1;
    t <<= 1;
  }

  while(m > 0) {
    m >>= 2;
    for(size_t j = 0; j < m; j++) {
      const uint64_t   m1  = 2 * (m + j);
      const uint64_t   m2  = 2 * m1;
      const uint64_t   k   = 4 * t * j;
      const loop_ctx_t ctx = {.w1     = {w[m1], w_con[m1]},
                              .r1     = {w[m2 + 0], w[m2 + 2]},
                              .r2     = {w[m2 + 1], w[m2 + 3]},
                              .r1_con = {w_con[m2 + 0], w_con[m2 + 2]},
                              .r2_con = {w_con[m2 + 1], w_con[m2 + 3]},
                              .zero   = {0},
                              .neg_q  = (-1 * q) & VMSL_WORD_SIZE_MASK,
                              .q2     = 2 * q,
                              .q4     = 4 * q,
                              .q2_vec = {2 * q, 2 * q},
                              .q4_vec = {4 * q, 4 * q}};

      if(t == 1) {
        for(size_t i = k; i < k + t; i++) {
          single_inv_butterfly(a, q, &ctx, i, t);
        }
      } else {
        for(size_t i = k; i < k + t; i += 2) {
          double_inv_butterfly(a, &ctx, i, t >> 1);
        }
      }
    }

    t <<= 2;
  }

  for(size_t i = 0; i < N; i++) {
    // At the last iteration, multiply by n^-1 mod q
    a[i] = fast_mul_mod_q(n_inv, a[i], q);
  }
}

/******************************
       Double input
******************************/
void fwd_ntt_radix4_intrinsic_lazy_dbl(uint64_t       a1[],
                                       uint64_t       a2[],
                                       const uint64_t N,
                                       const uint64_t q,
                                       const uint64_t w[],
                                       const uint64_t w_con[])
{
  const size_t bound_r4 = HAS_AN_EVEN_POWER(N) ? N : (N >> 1);

  for(size_t m = 1, t = (N >> 2); m < bound_r4; m <<= 2, t >>= 2) {
    for(size_t j = 0; j < m; j++) {
      const size_t     m1  = 2 * (m + j);
      const size_t     m2  = 2 * m1;
      const size_t     k   = 4 * t * j;
      const loop_ctx_t ctx = {.w1     = {w[m1], w_con[m1]},
                              .r1     = {w[m2 + 0], w[m2 + 1]},
                              .r2     = {w[m2 + 2], w[m2 + 3]},
                              .r1_con = {w_con[m2 + 0], w_con[m2 + 1]},
                              .r2_con = {w_con[m2 + 2], w_con[m2 + 3]},
                              .zero   = {0},
                              .neg_q  = (-1 * q) & VMSL_WORD_SIZE_MASK,
                              .q2     = 2 * q,
                              .q4     = 4 * q,
                              .q2_vec = {2 * q, 2 * q},
                              .q4_vec = {4 * q, 4 * q}};

      if(t == 1) {
        for(size_t i = k; i < k + t; i++) {
          single_fwd_butterfly(a1, q, &ctx, i, t);
          single_fwd_butterfly(a2, q, &ctx, i, t);
        }
      } else {
        for(size_t i = k; i < k + t; i += 2) {
          double_fwd_butterfly(a1, &ctx, i, t);
          double_fwd_butterfly(a2, &ctx, i, t);
        }
      }
    }
  }

  // Check whether N=2^m where m is odd.
  if(HAS_AN_EVEN_POWER(N)) {
    return;
  }

  for(size_t i = 0; i < N; i += 2) {
    const mul_op_t w1 = {w[i + N], w_con[i + N]};
    a1[i]             = reduce_8q_to_4q(a1[i], q);
    a2[i]             = reduce_8q_to_4q(a2[i], q);

    harvey_fwd_butterfly(&a1[i], &a1[i + 1], w1, q);
    harvey_fwd_butterfly(&a2[i], &a2[i + 1], w1, q);
  }
}
