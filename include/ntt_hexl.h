// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "defs.h"

EXTERNC_BEGIN

#ifdef AVX512_IFMA_SUPPORT

// Internal function of Intel HEXL under the license of Intel HEXL
void ForwardTransformToBitReverseAVX512(
  uint64_t *      operand,
  uint64_t        degree,
  uint64_t        mod,
  const uint64_t *root_of_unity_powers,
  const uint64_t *precon_root_of_unity_powers,
  uint64_t        input_mod_factor,
  uint64_t        output_mod_factor,
  uint64_t        recursion_depth,
  uint64_t        recursion_half);

static inline void fwd_ntt_radix2_hexl_lazy(uint64_t       a[],
                                            const uint64_t N,
                                            const uint64_t q,
                                            const uint64_t w[],
                                            const uint64_t w_con[])
{
  ForwardTransformToBitReverseAVX512(a, N, q, w, w_con, 2, 2, 0, 0);
}

static inline void fwd_ntt_radix2_hexl(uint64_t       a[],
                                       const uint64_t N,
                                       const uint64_t q,
                                       const uint64_t w[],
                                       const uint64_t w_con[])
{
  ForwardTransformToBitReverseAVX512(a, N, q, w, w_con, 2, 1, 0, 0);
}

#endif

EXTERNC_END
