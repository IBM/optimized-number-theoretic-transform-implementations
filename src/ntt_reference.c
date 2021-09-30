// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ntt_reference.h"
#include "fast_mul_operators.h"

/******************************
       Single input
******************************/

void fwd_ntt_ref_harvey_lazy(uint64_t       a[],
                             const uint64_t N,
                             const uint64_t q,
                             const uint64_t w[],
                             const uint64_t w_con[])
{
  size_t t = N >> 1;

  for(size_t m = 1; m < N; m <<= 1, t >>= 1) {
    size_t k = 0;
    for(size_t i = 0; i < m; i++) {
      const mul_op_t w1 = {w[m + i], w_con[m + i]};

      LOOP_UNROLL_4
      for(size_t j = k; j < k + t; j++) {
        harvey_fwd_butterfly(&a[j], &a[j + t], w1, q);
      }
      k = k + (2 * t);
    }
  }
}

void inv_ntt_ref_harvey(uint64_t       a[],
                        const uint64_t N,
                        const uint64_t q,
                        const mul_op_t n_inv,
                        const uint64_t word_size,
                        const uint64_t w[],
                        const uint64_t w_con[])
{
  uint64_t t = 1;

  for(size_t m = N >> 1; m > 1; m >>= 1, t <<= 1) {
    size_t k = 0;
    for(size_t i = 0; i < m; i++) {
      const mul_op_t w1 = {w[m + i], w_con[m + i]};

      for(size_t j = k; j < k + t; j++) {
        harvey_bkw_butterfly(&a[j], &a[j + t], w1, q);
      }
      k = k + (2 * t);
    }
  }

  // Final round - the harvey_bkw_butterfly, where the output is multiplies by
  // n_inv. Here m=1, k=0, t=N/2.
  const __uint128_t tmp = fast_mul_mod_q2(n_inv, w[1], q);

  // We can speed up this integer devision by using barreto reduction.
  // However, as it happens only once we keep the code simple.
  const mul_op_t w1 = {tmp, (tmp << word_size) / q};

  for(size_t j = 0; j < t; j++) {
    harvey_bkw_butterfly_final(&a[j], &a[j + t], w1, n_inv, q);
  }
}

/******************************
       Double input
******************************/
void fwd_ntt_ref_harvey_lazy_dbl(uint64_t       a1[],
                                 uint64_t       a2[],
                                 const uint64_t N,
                                 const uint64_t q,
                                 const uint64_t w[],
                                 const uint64_t w_con[])
{
  uint64_t t = N >> 1;

  for(size_t m = 1; m < N; m <<= 1, t >>= 1) {
    size_t k = 0;
    for(size_t i = 0; i < m; i++) {
      const mul_op_t w1 = {w[m + i], w_con[m + i]};
      for(size_t j = k; j < k + t; j++) {
        harvey_fwd_butterfly(&a1[j], &a1[j + t], w1, q);
        harvey_fwd_butterfly(&a2[j], &a2[j + t], w1, q);
      }
      k = k + (2 * t);
    }
  }
}
