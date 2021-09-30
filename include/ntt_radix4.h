// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "fast_mul_operators.h"

EXTERNC_BEGIN

void fwd_ntt_radix4_lazy(uint64_t       a[],
                         uint64_t       N,
                         uint64_t       q,
                         const uint64_t w[],
                         const uint64_t w_con[]);

static inline void fwd_ntt_radix4(uint64_t       a[],
                                  const uint64_t N,
                                  const uint64_t q,
                                  const uint64_t w[],
                                  const uint64_t w_con[])
{
  fwd_ntt_radix4_lazy(a, N, q, w, w_con);

  // Final reduction
  for(size_t i = 0; i < N; i++) {
    a[i] = reduce_8q_to_q(a[i], q);
  }
}

void inv_ntt_radix4(uint64_t       a[],
                    uint64_t       N,
                    uint64_t       q,
                    mul_op_t       n_inv,
                    const uint64_t w[],
                    const uint64_t w_con[]);

EXTERNC_END
