// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <string.h>

#include "measurements.h"
#include "ntt_radix4.h"
#include "ntt_radix4x4.h"
#include "ntt_reference.h"
#include "ntt_seal.h"
#include "tests.h"
#include "utils.h"

#ifdef S390X
#  include "ntt_radix4_s390x_vef.h"
#endif

#ifdef AVX512_IFMA_SUPPORT
#  include "ntt_avx512_ifma.h"
#  include "ntt_hexl.h"
#endif

void report_test_fwd_perf_headers(void)
{
  printf("                     |            fwd                                  "
         "                  |          fwd-lazy\n");
  printf("-----------------------------------------------------------------------"
         "---------------------");
  printf("--------------------------------------\n");
  printf("  N                q");
  printf("  rad2-ref");
  printf(" rad2-SEAL");
  printf("      rad4");
  printf("    rad4x4");
#ifdef S390X
  printf(" rad4-vmsl");
#elif AVX512_IFMA_SUPPORT
  printf(" rad2-hexl");
  printf(" rad2-ifma");
  printf(" rad2-ifma2");
  printf(" r4r2-ifma");
  printf(" r216-ifma");
#endif
  printf("  rad2-dbl");
#ifdef S390X
  printf(" rad4v-dbl");
#endif

  printf("  rad2-ref");
  printf(" rad2-SEAL");
  printf("      rad4");
#ifdef S390X
  printf(" rad4-vmsl");
#endif
  printf("\n");
}

static inline void test_fwd_perf(const test_case_t *t,
                                 uint64_t *         a,
                                 uint64_t *         b,
                                 const uint64_t *   a_cpy)
{
  const uint64_t q = t->q;
  const uint64_t n = t->n;

  printf("%3.0lu 0x%14.0lx ", t->m, t->q);

  MEASURE(fwd_ntt_ref_harvey(a, n, q, t->w_powers.ptr, t->w_powers_con.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_seal(a, n, q, t->w_powers.ptr, t->w_powers_con.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_radix4(a, n, q, t->w_powers_r4.ptr, t->w_powers_con_r4.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_radix4x4(a, n, q, t->w_powers_r4.ptr, t->w_powers_con_r4.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

#ifdef S390X
  MEASURE(fwd_ntt_radix4_intrinsic(a, n, q, t->w_powers_r4.ptr,
                                   t->w_powers_con_r4_vmsl.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));
#elif AVX512_IFMA_SUPPORT
  MEASURE(fwd_ntt_radix2_hexl(a, t->n, t->q, t->w_powers_hexl.ptr,
                              t->w_powers_con_hexl.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_radix4_avx512_ifma(a, t->n, t->q,
                                     t->w_powers_r4_avx512_ifma.ptr,
                                     t->w_powers_con_r4_avx512_ifma.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_radix4_avx512_ifma_unordered(
    a, t->n, t->q, t->w_powers_r4_avx512_ifma_unordered.ptr,
    t->w_powers_con_r4_avx512_ifma_unordered.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_r4r2_avx512_ifma(a, t->n, t->q,
                                   t->w_powers_r4r2_avx512_ifma.ptr,
                                   t->w_powers_con_r4r2_avx512_ifma.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_r2_16_avx512_ifma(a, t->n, t->q,
                                    t->w_powers_r2_16_avx512_ifma.ptr,
                                    t->w_powers_con_r2_16_avx512_ifma.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));
#endif

  MEASURE(fwd_ntt_ref_harvey_dbl(a, b, t->n, t->q, t->w_powers.ptr,
                                 t->w_powers_con.ptr));
  memcpy(a, a_cpy, n * sizeof(uint64_t));
  memcpy(b, a_cpy, n * sizeof(uint64_t));

#ifdef S390X
  MEASURE(fwd_ntt_radix4_intrinsic_dbl(a, b, t->n, t->q, t->w_powers_r4.ptr,
                                       t->w_powers_con_r4_vmsl.ptr));
#endif

  memcpy(a, a_cpy, n * sizeof(uint64_t));
  memcpy(b, a_cpy, n * sizeof(uint64_t));

  MEASURE(
    fwd_ntt_ref_harvey_lazy(a, n, q, t->w_powers.ptr, t->w_powers_con.ptr););
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(fwd_ntt_seal_lazy(a, n, q, t->w_powers.ptr, t->w_powers_con.ptr););
  memcpy(a, a_cpy, n * sizeof(uint64_t));

  MEASURE(
    fwd_ntt_radix4_lazy(a, n, q, t->w_powers_r4.ptr, t->w_powers_con_r4.ptr););
  memcpy(a, a_cpy, n * sizeof(uint64_t));

#ifdef S390X
  MEASURE(fwd_ntt_radix4_intrinsic_lazy(a, n, q, t->w_powers_r4.ptr,
                                        t->w_powers_con_r4_vmsl.ptr));
#endif

  printf("\n");
}

void test_aligned_fwd_perf(const test_case_t *t)
{
  const uint64_t n = t->n;
  const uint64_t q = t->q;

  // We use a_cpy to reset a after every NTT call.
  // This is especially important when dealing with the lazy evaluation functions
  // To avoid overflowing and therefore slowdowns of VMSL.
  ALIGN(64) uint64_t a[n];
  ALIGN(64) uint64_t b[n];
  ALIGN(64) uint64_t a_cpy[n];
  random_buf(a, n, q);
  memcpy(a_cpy, a, sizeof(a));
  memcpy(b, a, sizeof(a));

  test_fwd_perf(t, a, b, a_cpy);
}

void test_unaligned_fwd_perf(const test_case_t *t)
{
  const uint64_t n = t->n;
  const uint64_t q = t->q;

  // We use a_cpy to reset a after every NTT call.
  // This is especially important when dealing with the lazy evaluation functions
  // To avoid overflowing and therefore slowdowns of VMSL.
  unaligned64_ptr_t a;
  unaligned64_ptr_t b;
  unaligned64_ptr_t a_cpy;
  allocate_unaligned_array(&a, n);
  allocate_unaligned_array(&b, n);
  allocate_unaligned_array(&a_cpy, n);

  random_buf(a.ptr, n, q);
  memset(a_cpy.ptr, 0, n * sizeof(uint64_t));
  memcpy(a_cpy.ptr, a.ptr, n * sizeof(uint64_t));
  memset(b.ptr, 0, n * sizeof(uint64_t));
  memcpy(b.ptr, a.ptr, n * sizeof(uint64_t));

  test_fwd_perf(t, a.ptr, b.ptr, a_cpy.ptr);

  free_unaligned_array(&a);
  free_unaligned_array(&b);
  free_unaligned_array(&a_cpy);
}

void report_test_inv_perf_headers(void)
{
  printf("                     |            inv\n");
  printf("------------------------------------------------------\n");
  printf("  N                q");

  printf("  rad2-ref");
  printf(" rad2-SEAL");
  printf("      rad4");

#ifdef S390X
  printf(" rad4-vmsl");
#endif

  printf("\n");
}

void test_inv_perf(const test_case_t *t)
{
  const uint64_t n = t->n;
  const uint64_t q = t->q;

  printf("%3.0lu 0x%14.0lx ", t->m, t->q);

  // We use a_cpy to reset a after every NTT call.
  // This is especially important when dealing with the lazy evaluation functions
  // To avoid overflowing and therefore slowdowns of VMSL.
  uint64_t a[n];
  uint64_t a_cpy[n];
  random_buf(a, n, q);
  memcpy(a_cpy, a, sizeof(a));

  MEASURE(inv_ntt_ref_harvey(a, n, q, t->n_inv, WORD_SIZE, t->w_inv_powers.ptr,
                             t->w_inv_powers_con.ptr));
  memcpy(a, a_cpy, sizeof(a));

  MEASURE(inv_ntt_seal(a, t->n, t->q, t->n_inv.op, t->n_inv.con,
                       t->w_inv_powers.ptr, t->w_inv_powers_con.ptr));
  memcpy(a, a_cpy, sizeof(a));

  MEASURE(inv_ntt_radix4(a, n, q, t->n_inv, t->w_inv_powers_r4.ptr,
                         t->w_inv_powers_con_r4.ptr));
  memcpy(a, a_cpy, sizeof(a));

#ifdef S390X
  MEASURE(inv_ntt_radix4_intrinsic(a, n, q, t->n_inv_vmsl, t->w_inv_powers_r4.ptr,
                                   t->w_inv_powers_con_r4_vmsl.ptr));
#endif

  printf("\n");
}

void test_fwd_single_case(const test_case_t *t, const func_num_t func_num)
{
  const uint64_t n = t->n;
  const uint64_t q = t->q;

  // We use a_cpy to reset a after every NTT call.
  // This is especially important when dealing with the lazy evaluation functions
  // To avoid overflowing and therefore slowdowns of VMSL.
  uint64_t a[n];
  uint64_t a_cpy[n];
  random_buf(a, n, q);
  memcpy(a_cpy, a, sizeof(a));

  switch(func_num) {
    case FWD_REF:
      MEASURE(fwd_ntt_ref_harvey(a, n, q, t->w_powers.ptr, t->w_powers_con.ptr));
      break;
    case FWD_SEAL:
      MEASURE(fwd_ntt_seal(a, n, q, t->w_powers.ptr, t->w_powers_con.ptr));
      break;
    case FWD_R4:
      MEASURE(
        fwd_ntt_radix4(a, n, q, t->w_powers_r4.ptr, t->w_powers_con_r4.ptr));
      break;
    case FWD_R4x4:
      MEASURE(
        fwd_ntt_radix4x4(a, n, q, t->w_powers_r4.ptr, t->w_powers_con_r4.ptr));
      break;
    case FWD_R4_VMSL:
#ifdef S390X
      MEASURE(fwd_ntt_radix4_intrinsic(a, n, q, t->w_powers_r4.ptr,
                                       t->w_powers_con_r4_vmsl.ptr));
      break;
#elif AVX512_IFMA_SUPPORT
    case FWD_R4_HEXL:
      MEASURE(fwd_ntt_radix2_hexl(a, t->n, t->q, t->w_powers_hexl.ptr,
                                  t->w_powers_con_hexl.ptr));
      break;
    case FWD_R4_AVX512_IFMA:
      MEASURE(fwd_ntt_radix4_avx512_ifma(a, t->n, t->q,
                                         t->w_powers_r4_avx512_ifma.ptr,
                                         t->w_powers_con_r4_avx512_ifma.ptr));
      break;
    case FWD_R4_AVX512_IFMA_UNORDERED:
      MEASURE(fwd_ntt_radix4_avx512_ifma_unordered(
        a, t->n, t->q, t->w_powers_r4_avx512_ifma_unordered.ptr,
        t->w_powers_con_r4_avx512_ifma_unordered.ptr));
      break;
    case FWD_R4R2_AVX512_IFMA:
      MEASURE(fwd_ntt_r4r2_avx512_ifma(a, t->n, t->q,
                                       t->w_powers_r4r2_avx512_ifma.ptr,
                                       t->w_powers_con_r4r2_avx512_ifma.ptr));
      break;
    case FWD_R2_R16_AVX512_IFMA:
      MEASURE(fwd_ntt_r2_16_avx512_ifma(a, t->n, t->q,
                                        t->w_powers_r2_16_avx512_ifma.ptr,
                                        t->w_powers_con_r2_16_avx512_ifma.ptr));
      break;
#endif
    default: break;
  }
}
