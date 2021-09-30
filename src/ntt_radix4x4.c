// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ntt_radix4x4.h"
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

static inline uint64_t get_iter_reminder(const uint64_t N)
{
  if(HAS_AN_REM1_POWER(N)) {
    return 1;
  }
  if(HAS_AN_REM2_POWER(N)) {
    return 2;
  }
  if(HAS_AN_REM3_POWER(N)) {
    return 3;
  }
  return 0;
}

void fwd_ntt_radix4x4_lazy(uint64_t       a[],
                           const uint64_t N,
                           const uint64_t q,
                           const uint64_t w[],
                           const uint64_t w_con[])
{
  const uint64_t bound_r4 = N;
  uint64_t       m_rem    = get_iter_reminder(N);

  mul_op_t roots[5];
  mul_op_t roots4[4][5];
  size_t   t = N >> 2;

  for(size_t m = 1; m < (bound_r4 >> m_rem); m <<= 4) {
    for(size_t j = 0; j < m; j++) {
      const uint64_t k  = 4 * t * j;
      size_t         t2 = t >> 2;

      collect_roots(roots, w, w_con, m, j);
      for(size_t i = 0; i < 4; i++) {
        collect_roots(roots4[i], w, w_con, m << 2, 4 * j + i);
      }

      // Perform the 16-radix NTT in two steps of radix-4 NTT
      for(size_t i = k; i < k + t2; i++) {
        for(size_t l = i; l < i + t; l += t2) {
          radix4_fwd_butterfly(&a[l], &a[l + t], &a[l + 2 * t], &a[l + 3 * t],
                               roots, q);
        }
        size_t x = 0;
        for(size_t l = i; l < i + 4 * t; l += t, x++) {
          radix4_fwd_butterfly(&a[l], &a[l + t2], &a[l + 2 * t2], &a[l + 3 * t2],
                               roots4[x], q);
        }
      }
    }
    t >>= 4;
  }

  // Perform extra iterations if needed
  switch(m_rem) {
    case 1:
      // Perform extra radix-2 iteration.
      for(size_t i = 0; i < N; i += 2) {
        const mul_op_t w1 = {w[N + i], w_con[N + i]};
        a[i]              = reduce_8q_to_4q(a[i], q);

        harvey_fwd_butterfly(&a[i], &a[i + 1], w1, q);
      }
      return;
    case 3:
      // Perform extra radix-2 and then radix-4 iteration.
      t              = 4;
      const size_t m = N >> 3;
      for(size_t i = 0; i < m; i++) {
        const size_t   k  = 2 * t * i;
        const mul_op_t w1 = {w[2 * (m + i)], w_con[2 * (m + i)]};
        a[i]              = reduce_8q_to_4q(a[i], q);

        for(size_t j = k; j < k + t; j++) {
          harvey_fwd_butterfly(&a[j], &a[j + t], w1, q);
        }
      }
      /* fall through */
    case 2:
      // Perform extra radix-4 iteration (for cases 2 and 3).
      for(size_t i = 0; i < N; i += 4) {
        collect_roots(roots, w, w_con, N >> 2, i >> 2);
        radix4_fwd_butterfly(&a[i], &a[i + 1], &a[i + 2], &a[i + 3], roots, q);
      }
      return;
    default: return;
  }
}
