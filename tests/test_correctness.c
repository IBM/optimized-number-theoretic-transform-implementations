// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <string.h>

#include "ntt_radix4.h"
#include "ntt_radix4x4.h"
#include "ntt_reference.h"
#include "ntt_seal.h"
#include "pre_compute.h"
#include "test_cases.h"
#include "utils.h"

#ifdef S390X
#  include "ntt_radix4_s390x_vef.h"
#endif

#ifdef AVX512_IFMA_SUPPORT
#  include "ntt_avx512_ifma.h"
#  include "ntt_hexl.h"
#endif

static inline int test_radix2_scalar(const test_case_t *t, uint64_t a_orig[])
{
  uint64_t a[t->n];
  memcpy(a, a_orig, sizeof(a));

  printf("Running fwd_ntt_ref_harvey\n");
  fwd_ntt_ref_harvey(a, t->n, t->q, t->w_powers, t->w_powers_con);

  printf("Running inv_ntt_ref_harvey\n");
  inv_ntt_ref_harvey(a, t->n, t->q, t->n_inv, WORD_SIZE, t->w_inv_powers,
                     t->w_inv_powers_con);

  GUARD_MSG(memcmp(a_orig, a, sizeof(a)), "Bad results after radix-2 inv\n");

  return SUCCESS;
}

static inline int test_radix2_scalar_dbl(const test_case_t *t,
                                         uint64_t           a_orig[],
                                         uint64_t           b_orig[],
                                         uint64_t           a_ntt[])
{
  uint64_t a[t->n];
  uint64_t b[t->n];
  memcpy(a, a_orig, sizeof(a));
  memcpy(b, b_orig, sizeof(a));

  printf("Running fwd_ntt_ref_harvey_dbl\n");
  fwd_ntt_ref_harvey_dbl(a, b, t->n, t->q, t->w_powers, t->w_powers_con);

  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after radix-2 scalar double for a\n");
  GUARD_MSG(memcmp(a_ntt, b, sizeof(a)),
            "Bad results after radix-2 scalar double for b\n");

  return SUCCESS;
}

static inline int
test_radix2_scalar_seal(const test_case_t *t, uint64_t a_orig[], uint64_t a_ntt[])
{
  uint64_t a[t->n];
  memcpy(a, a_orig, sizeof(a));

  printf("Running fwd_ntt_seal\n");
  fwd_ntt_seal(a, t->n, t->q, t->w_powers, t->w_powers_con);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after radix-2 SEAL fwd implementation\n");

  printf("Running inv_ntt_seal\n");
  inv_ntt_seal(a, t->n, t->q, t->n_inv.op, t->n_inv.con, t->w_inv_powers,
               t->w_inv_powers_con);
  GUARD_MSG(memcmp(a_orig, a, sizeof(a)),
            "Bad results after radix-2 SEAL inv implementation\n");

  return SUCCESS;
}

static inline int
test_radix4_scalar(const test_case_t *t, uint64_t a_orig[], uint64_t a_ntt[])
{
  uint64_t a[t->n];
  memcpy(a, a_orig, sizeof(a));

  printf("Running fwd_ntt_radix4\n");
  fwd_ntt_radix4(a, t->n, t->q, t->w_powers_r4, t->w_powers_con_r4);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)), "Bad results after radix-4 fwd\n");

  printf("Running inv_ntt_radix4\n");
  inv_ntt_radix4(a, t->n, t->q, t->n_inv, t->w_inv_powers_r4,
                 t->w_inv_powers_con_r4);

  GUARD_MSG(memcmp(a_orig, a, sizeof(a)), "Bad results after radix-4 inv\n");

  return SUCCESS;
}

static inline int
test_radix4x4_scalar(const test_case_t *t, uint64_t a_orig[], uint64_t a_ntt[])
{
  uint64_t a[t->n];
  memcpy(a, a_orig, sizeof(a));

  printf("Running fwd_ntt_radix4x4\n");
  fwd_ntt_radix4x4(a, t->n, t->q, t->w_powers_r4, t->w_powers_con_r4);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)), "Bad results after radix-4x4 fwd\n");
  /*
    printf("Running inv_ntt_radix4\n");
    inv_ntt_radix4(a, t->n, t->q, t->n_inv, t->w_inv_powers_r4,
                   t->w_inv_powers_con_r4);

    GUARD_MSG(memcmp(a_orig, a, sizeof(a)), "Bad results after radix-4 inv\n");
  */
  return SUCCESS;
}

#ifdef S390X
static inline int
test_radix4_intrinsic(const test_case_t *t, uint64_t a_orig[], uint64_t a_ntt[])
{
  uint64_t a[t->n];
  memcpy(a, a_orig, sizeof(a));

  printf("Running fwd_ntt_radix4_intrinsic\n");
  fwd_ntt_radix4_intrinsic(a, t->n, t->q, t->w_powers_r4,
                           t->w_powers_con_r4_vmsl);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after radix-4 with intrinsic fwd\n");

  printf("Running inv_ntt_radix4_intrinsic\n");
  inv_ntt_radix4_intrinsic(a, t->n, t->q, t->n_inv_vmsl, t->w_inv_powers_r4,
                           t->w_inv_powers_con_r4_vmsl);

  GUARD_MSG(memcmp(a_orig, a, sizeof(a)),
            "Bad results after radix-4 inv with intrinsic\n");

  return SUCCESS;
}

static inline int test_radix4_intrinsic_dbl(const test_case_t *t,
                                            uint64_t           a_orig[],
                                            uint64_t           b_orig[],
                                            uint64_t           a_ntt[])
{
  uint64_t a[t->n];
  uint64_t b[t->n];
  memcpy(a, a_orig, sizeof(a));
  memcpy(b, b_orig, sizeof(a));

  printf("Running fwd_ntt_ref_harvey_dbl\n");
  fwd_ntt_radix4_intrinsic_dbl(a, b, t->n, t->q, t->w_powers_r4,
                               t->w_powers_con_r4_vmsl);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after radix-2 scalar double for a\n");
  GUARD_MSG(memcmp(a_ntt, b, sizeof(b)),
            "Bad results after radix-2 scalar double for b\n");

  return SUCCESS;
}
#endif

#ifdef AVX512_IFMA_SUPPORT
static inline int
test_radix2_hexl(const test_case_t *t, uint64_t a_orig[], uint64_t a_ntt[])
{
  // We can't test AVX512-IFMA for q>2^49 so we always success.
  if(t->q & AVX512_IFMA_MAX_MODULUS_MASK) {
    return SUCCESS;
  }

  uint64_t a[t->n];
  memcpy(a, a_orig, sizeof(a));

  printf("Running fwd_ntt_radix2_hexl\n");
  fwd_ntt_radix2_hexl(a, t->n, t->q, t->w_powers_hexl, t->w_powers_con_hexl);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after HEXL radix-2 with AVX512-IFMA intrinsic fwd\n");

  return SUCCESS;
}

static inline void fix_a_order(uint64_t *a, uint64_t n)
{
  const __m512i idx = _mm512_setr_epi64(0, 4, 8, 12, 16, 20, 24, 28);

  for(size_t i = 0; i < n; i += (4 * 8)) {
    __m512i X = LOAD(&a[i + 0]);
    __m512i Y = LOAD(&a[i + 8]);
    __m512i Z = LOAD(&a[i + 16]);
    __m512i T = LOAD(&a[i + 24]);

    SCATTER(&a[i + 0], idx, X, 8);
    SCATTER(&a[i + 1], idx, Y, 8);
    SCATTER(&a[i + 2], idx, Z, 8);
    SCATTER(&a[i + 3], idx, T, 8);

    X = LOAD(&a[i + 0]);
    Y = LOAD(&a[i + 8]);
    Z = LOAD(&a[i + 16]);
    T = LOAD(&a[i + 24]);

    const __m512i X1 = SHUF(X, Y, 0x44);
    const __m512i Z1 = SHUF(X, Y, 0xee);
    const __m512i Y1 = SHUF(Z, T, 0x44);
    const __m512i T1 = SHUF(Z, T, 0xee);

    STORE(&a[i + 0], X1);
    STORE(&a[i + 8], Y1);
    STORE(&a[i + 16], Z1);
    STORE(&a[i + 24], T1);
  }
}

static inline int
test_radix4_avx512_ifma(const test_case_t *t, uint64_t a_orig[], uint64_t a_ntt[])
{
  // We can't test AVX512-IFMA for q>2^49 so we always success.
  if(t->q & AVX512_IFMA_MAX_MODULUS_MASK) {
    return SUCCESS;
  }

  uint64_t a[t->n];
  memcpy(a, a_orig, sizeof(a));

  printf("Running fwd_ntt_radix4_avx512_ifma\n");
  fwd_ntt_radix4_avx512_ifma(a, t->n, t->q, t->w_powers_r4_avx512_ifma,
                             t->w_powers_con_r4_avx512_ifma);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after radix-4 with AVX512-IFMA intrinsic fwd\n");

  memcpy(a, a_orig, sizeof(a));
  printf("Running fwd_ntt_radix4_avx512_ifma_unordered\n");
  fwd_ntt_radix4_avx512_ifma_unordered(a, t->n, t->q,
                                       t->w_powers_r4_avx512_ifma_unordered,
                                       t->w_powers_con_r4_avx512_ifma_unordered);
  fix_a_order(a, t->n);
  GUARD_MSG(
    memcmp(a_ntt, a, sizeof(a)),
    "Bad results after radix-4 with AVX512-IFMA intrinsic unordered fwd\n");

  memcpy(a, a_orig, sizeof(a));
  printf("Running fwd_ntt_r4r2_avx512_ifma\n");
  fwd_ntt_r4r2_avx512_ifma(a, t->n, t->q, t->w_powers_r4r2_avx512_ifma,
                           t->w_powers_con_r4r2_avx512_ifma);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after r4r2 with AVX512-IFMA intrinsic fwd\n");

  memcpy(a, a_orig, sizeof(a));
  printf("Running fwd_ntt_r2_16_avx512_ifma\n");
  fwd_ntt_r2_16_avx512_ifma(a, t->n, t->q, t->w_powers_r2_16_avx512_ifma,
                            t->w_powers_con_r2_16_avx512_ifma);
  GUARD_MSG(memcmp(a_ntt, a, sizeof(a)),
            "Bad results after r2_16 with AVX512-IFMA intrinsic fwd\n");

  return SUCCESS;
}
#endif

int test_correctness(const test_case_t *t)
{
  // Prepare input
  uint64_t a[t->n];
  uint64_t b[t->n];
  uint64_t a_ntt[t->n];
  uint64_t a_cpy[t->n];
  random_buf(a, t->n, t->q);
  memcpy(a_cpy, a, sizeof(a));
  memcpy(b, a, sizeof(a));

  // Prepare a_ntt = NTT(a)
  fwd_ntt_ref_harvey(a_cpy, t->n, t->q, t->w_powers, t->w_powers_con);
  memcpy(a_ntt, a_cpy, sizeof(a_cpy));

  GUARD(test_radix2_scalar(t, a));
  GUARD(test_radix2_scalar_dbl(t, a, b, a_ntt));
  GUARD(test_radix2_scalar_seal(t, a, a_ntt))
  GUARD(test_radix4_scalar(t, a, a_ntt))
  GUARD(test_radix4x4_scalar(t, a, a_ntt))
#ifdef S390X
  GUARD(test_radix4_intrinsic(t, a, a_ntt))
  GUARD(test_radix4_intrinsic_dbl(t, a, b, a_ntt))
#elif AVX512_IFMA_SUPPORT
  GUARD(test_radix2_hexl(t, a, a_ntt))
  GUARD(test_radix4_avx512_ifma(t, a, a_ntt))
#endif

  return SUCCESS;
}
