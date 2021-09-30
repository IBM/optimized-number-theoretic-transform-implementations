// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "fast_mul_operators.h"

EXTERNC_BEGIN

/******************************
       Single input
******************************/
void fwd_ntt_ref_harvey_lazy(uint64_t       a[],
                             uint64_t       N,
                             uint64_t       q,
                             const uint64_t w[],
                             const uint64_t w_con[]);

static inline void fwd_ntt_ref_harvey(uint64_t       a[],
                                      const uint64_t N,
                                      const uint64_t q,
                                      const uint64_t w[],
                                      const uint64_t w_con[])
{
  fwd_ntt_ref_harvey_lazy(a, N, q, w, w_con);

  // Final reduction
  for(size_t i = 0; i < N; i++) {
    a[i] = reduce_4q_to_q(a[i], q);
  }
}

void inv_ntt_ref_harvey(uint64_t       a[],
                        uint64_t       N,
                        uint64_t       q,
                        mul_op_t       n_inv,
                        uint64_t       word_size,
                        const uint64_t w[],
                        const uint64_t w_con[]);

/******************************
       Double input
******************************/
void fwd_ntt_ref_harvey_lazy_dbl(uint64_t       a1[],
                                 uint64_t       a2[],
                                 uint64_t       N,
                                 uint64_t       q,
                                 const uint64_t w[],
                                 const uint64_t w_con[]);

static inline void fwd_ntt_ref_harvey_dbl(uint64_t       a1[],
                                          uint64_t       a2[],
                                          const uint64_t N,
                                          const uint64_t q,
                                          const uint64_t w[],
                                          const uint64_t w_con[])
{
  fwd_ntt_ref_harvey_lazy_dbl(a1, a2, N, q, w, w_con);

  // Final reduction
  for(size_t i = 0; i < N; i++) {
    a1[i] = reduce_4q_to_q(a1[i], q);
    a2[i] = reduce_4q_to_q(a2[i], q);
  }
}

EXTERNC_END
