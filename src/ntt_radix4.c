// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ntt_radix4.h"
#include "fast_mul_operators.h"

static inline void collect_roots(mul_op_t       w1[5],
                                 const uint64_t w[],
                                 const uint64_t w_con[],
                                 const size_t   m,
                                 const size_t   j)
{
  const uint64_t m1 = 2 * (m + j);
  w1[0].op          = w[m1];
  w1[1].op          = w[2 * m1];
  w1[2].op          = w[2 * m1 + 1];
  w1[3].op          = w[2 * m1 + 2];
  w1[4].op          = w[2 * m1 + 3];

  w1[0].con = w_con[m1];
  w1[1].con = w_con[2 * m1];
  w1[2].con = w_con[2 * m1 + 1];
  w1[3].con = w_con[2 * m1 + 2];
  w1[4].con = w_con[2 * m1 + 3];
}

void fwd_ntt_radix4_lazy(uint64_t       a[],
                         const uint64_t N,
                         const uint64_t q,
                         const uint64_t w[],
                         const uint64_t w_con[])
{
  const uint64_t bound_r4 = HAS_AN_EVEN_POWER(N) ? N : (N >> 1);
  mul_op_t       roots[5];
  size_t         t = N >> 2;

  for(size_t m = 1; m < bound_r4; m <<= 2) {
    for(size_t j = 0; j < m; j++) {
      const uint64_t k = 4 * t * j;

      collect_roots(roots, w, w_con, m, j);
      for(size_t i = k; i < k + t; i++) {
        radix4_fwd_butterfly(&a[i], &a[i + t], &a[i + 2 * t], &a[i + 3 * t],
                             roots, q);
      }
    }
    t >>= 2;
  }

  // Check whether N=2^m where m is odd.
  // If not perform extra radix-2 iteration.
  if(HAS_AN_EVEN_POWER(N)) {
    return;
  }

  for(size_t i = 0; i < N; i += 2) {
    const mul_op_t w1 = {w[N + i], w_con[N + i]};
    a[i]              = reduce_8q_to_4q(a[i], q);

    harvey_fwd_butterfly(&a[i], &a[i + 1], w1, q);
  }
}

void inv_ntt_radix4(uint64_t       a[],
                    const uint64_t N,
                    const uint64_t q,
                    const mul_op_t n_inv,
                    const uint64_t w[],
                    const uint64_t w_con[])
{
  uint64_t t = 1;
  uint64_t m = N;
  mul_op_t roots[5];

  // 1. Check whether N=2^m where m is even.
  // If yes, reduce all values modulo 2q, this also can be done outside of this
  // function. Otherwise, perform one radix-2 iteration.
  if(HAS_AN_EVEN_POWER(N)) {
    for(size_t i = 0; i < N; i++) {
      a[i] = reduce_8q_to_2q(a[i], q);
    }

  } else {
    // Perform the first iteration as a radix-2 iteration.
    for(size_t i = 0; i < N; i += 2) {
      const mul_op_t w1 = {w[N + i], w_con[N + i]};

      a[i] = reduce_8q_to_4q(a[i], q);
      harvey_bkw_butterfly(&a[i], &a[i + 1], w1, q);
    }

    m >>= 1;
    t <<= 1;
  }

  // 2. Perform radix-4 NTT iterations.
  for(m >>= 2; m > 0; m >>= 2) {
    for(size_t j = 0; j < m; j++) {
      const uint64_t k = 4 * t * j;
      collect_roots(roots, w, w_con, m, j);

      for(size_t i = k; i < k + t; i++) {
        radix4_inv_butterfly(&a[i], &a[i + t], &a[i + 2 * t], &a[i + 3 * t],
                             roots, q);
      }
    }
    t <<= 2;
  }

  // 3. Normalize the results
  for(size_t i = 0; i < N; i++) {
    a[i] = fast_mul_mod_q(n_inv, a[i], q);
  }
}
