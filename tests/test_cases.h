// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdlib.h>

#include "fast_mul_operators.h"
#include "pre_compute.h"

EXTERNC_BEGIN

typedef struct test_case_s {
  // These parameters are predefined
  uint64_t m;
  uint64_t q;
  uint64_t w;
  uint64_t w_inv; // w^(-1) mod q
  mul_op_t n_inv; // 2^(-m) mod q

  // These parameters are dinamically computed based on the above values.
  uint64_t  n;
  uint64_t  qneg;
  uint64_t  q2;
  uint64_t  q4;
  uint64_t *w_powers;
  uint64_t *w_powers_con;
  uint64_t *w_inv_powers;
  uint64_t *w_inv_powers_con;

  // For radix-4 tests
  uint64_t *w_powers_r4;
  uint64_t *w_powers_con_r4;
  uint64_t *w_inv_powers_r4;
  uint64_t *w_inv_powers_con_r4;

#ifdef S390X
  // For radix-4 tests with VMSL (56-bits instead of 64-bits)
  uint64_t *w_powers_con_r4_vmsl;
  uint64_t *w_inv_powers_con_r4_vmsl;
  mul_op_t  n_inv_vmsl;
#endif

#ifdef AVX512_IFMA_SUPPORT
  // For radix-2 tests with AVX512-IFMA on X86-64 bit platofrms (52-bits)
  uint64_t *w_powers_hexl;
  uint64_t *w_powers_con_hexl;

  uint64_t *w_powers_r4_avx512_ifma;
  uint64_t *w_powers_con_r4_avx512_ifma;

  uint64_t *w_powers_r4_avx512_ifma_unordered;
  uint64_t *w_powers_con_r4_avx512_ifma_unordered;

  uint64_t *w_powers_r4r2_avx512_ifma;
  uint64_t *w_powers_con_r4r2_avx512_ifma;

  uint64_t *w_powers_r2_16_avx512_ifma;
  uint64_t *w_powers_con_r2_16_avx512_ifma;
#endif

} test_case_t;

/*
We used the following sagemath script to generate the values below.

params = [(7681, 8), (65537, 9),
          (65537, 10), (65537, 11),
          (65537, 12), (65537, 13),
          (65537, 14), (0xc0001, 14),
          (0xfff0001, 14), (0x1ffc8001, 14),
          (0x7ffe0001, 14), (0xfff88001, 14),
          (0x7fffffffe0001, 14), (0x80000001c0001, 14),
          (0x80000001c0001, 15),(0x7fffffffe0001, 16)
         ]

for (p, m) in params:
    n=2^m
    Zp = Integers(p)
    w_p = primitive_root(p)

    # Get a primitive 2n'th root of unity w
    eta = (p-1)//(2*n)
    w = Zp(w_p)^eta

    # Find minimum root
    new_w, min_w = w, w
    for i in range(2*n):
        min_w = min(min_w, new_w)
        new_w = new_w * w^2;
    print("p=",hex(p), " m=",m, " w=", min_w, " verify=", 1 == min_w^(2*n), "
w_inv=", min_w^-1, " n_inv=",Zp(n)^-1, sep='')
*/

// We use NOLINT in order to stop clang-tidy from reporting magic numbers
static test_case_t tests[] = {
  {.m = 8, .q = 0x1e01, .w = 62, .w_inv = 1115, .n_inv.op = 7651},      // NOLINT
  {.m = 9, .q = 0x10001, .w = 431, .w_inv = 55045, .n_inv.op = 65409},  // NOLINT
  {.m = 10, .q = 0x10001, .w = 33, .w_inv = 1986, .n_inv.op = 65473},   // NOLINT
  {.m = 11, .q = 0x10001, .w = 21, .w_inv = 49933, .n_inv.op = 65505},  // NOLINT
  {.m = 12, .q = 0x10001, .w = 13, .w_inv = 15124, .n_inv.op = 65521},  // NOLINT
  {.m = 13, .q = 0x10001, .w = 15, .w_inv = 30584, .n_inv.op = 65529},  // NOLINT
  {.m = 14, .q = 0x10001, .w = 9, .w_inv = 7282, .n_inv.op = 65533},    // NOLINT
  {.m = 14, .q = 0xc0001, .w = 9, .w_inv = 174763, .n_inv.op = 786385}, // NOLINT
  {.m        = 14,
   .q        = 0xfff0001,  // NOLINT
   .w        = 10360,      // NOLINT
   .w_inv    = 28987060,   // NOLINT
   .n_inv.op = 268353541}, // NOLINT
  {.m        = 14,
   .q        = 0x1ffc8001, // NOLINT
   .w        = 101907,
   .w_inv    = 42191135,   // NOLINT
   .n_inv.op = 536608783}, // NOLINT
  {.m        = 14,
   .q        = 0x7ffe0001, // NOLINT
   .w        = 320878,
   .w_inv    = 74168714,    // NOLINT
   .n_inv.op = 2147221513}, // NOLINT
  {.m        = 14,
   .q        = 0xfff88001, // NOLINT
   .w        = 263641,
   .w_inv    = 243522111,   // NOLINT
   .n_inv.op = 4294213663}, // NOLINT
  {.m        = 14,
   .q        = 0x7fffffffe0001, // NOLINT
   .w        = 83051296654,
   .w_inv    = 374947202223591,   // NOLINT
   .n_inv.op = 2251662374600713}, // NOLINT
  {.m        = 14,
   .q        = 0x80000001c0001, // NOLINT
   .w        = 72703961923,
   .w_inv    = 153477749218715,   // NOLINT
   .n_inv.op = 2251662376566673}, // NOLINT
  {.m        = 15,
   .q        = 0x80000001c0001, // NOLINT
   .w        = 82138512871,
   .w_inv    = 535648572761016,   // NOLINT
   .n_inv.op = 2251731096043465}, // NOLINT
  {.m        = 16,                // NOLINT
   .q        = 0x7fffffffe0001,   // NOLINT
   .w        = 29454831443,
   .w_inv    = 520731633805630,    // NOLINT
   .n_inv.op = 2251765453815811}}; // NOLINT

#define NUM_OF_TEST_CASES (sizeof(tests) / sizeof(test_case_t))

#define _ALLOCATE_U64_ARRAY(ptr, size)                                    \
  do {                                                                    \
    if(NULL == ((ptr) = (uint64_t *)malloc((size) * sizeof(uint64_t)))) { \
      printf("Allocation error");                                         \
      return 0;                                                           \
    }                                                                     \
  } while(0)

static inline int _init_test(test_case_t *t)
{
  // For brevity
  const uint64_t q     = t->q;
  const uint64_t w     = t->w;
  const uint64_t m     = t->m;
  const uint64_t w_inv = t->w_inv;
  const uint64_t n     = 1UL << t->m;

  t->n         = n;
  t->n_inv.con = calc_ninv_con(t->n_inv.op, q, WORD_SIZE);
  t->q2        = 2 * q;
  t->q4        = 4 * q;

  // Prepare radix-2 w-powers
  _ALLOCATE_U64_ARRAY(t->w_powers, n);
  calc_w(t->w_powers, w, n, q, m);

  _ALLOCATE_U64_ARRAY(t->w_powers_con, n);
  calc_w_con(t->w_powers_con, t->w_powers, n, q, WORD_SIZE);

  _ALLOCATE_U64_ARRAY(t->w_inv_powers, n);
  calc_w_inv(t->w_inv_powers, w_inv, n, q, m);

  _ALLOCATE_U64_ARRAY(t->w_inv_powers_con, n);
  calc_w_con(t->w_inv_powers_con, t->w_inv_powers, n, q, WORD_SIZE);

  // Expand the list of powers to support the radix-4 case.
  _ALLOCATE_U64_ARRAY(t->w_powers_r4, 2 * n);
  expand_w(t->w_powers_r4, t->w_powers, n, q);

  _ALLOCATE_U64_ARRAY(t->w_powers_con_r4, 2 * n);
  calc_w_con(t->w_powers_con_r4, t->w_powers_r4, 2 * n, q, WORD_SIZE);

  _ALLOCATE_U64_ARRAY(t->w_inv_powers_r4, 2 * n);
  expand_w(t->w_inv_powers_r4, t->w_inv_powers, n, q);

  _ALLOCATE_U64_ARRAY(t->w_inv_powers_con_r4, 2 * n);
  calc_w_con(t->w_inv_powers_con_r4, t->w_inv_powers_r4, 2 * n, q, WORD_SIZE);

#ifdef S390X
  t->n_inv_vmsl.con = calc_ninv_con(t->n_inv.op, q, VMSL_WORD_SIZE);
  t->n_inv_vmsl.op  = t->n_inv.op;

  // for radix-4 vmsl
  _ALLOCATE_U64_ARRAY(t->w_powers_con_r4_vmsl, 2 * n);
  calc_w_con(t->w_powers_con_r4_vmsl, t->w_powers_r4, 2 * t->n, t->q,
             VMSL_WORD_SIZE);

  _ALLOCATE_U64_ARRAY(t->w_inv_powers_con_r4_vmsl, 2 * n);
  calc_w_con(t->w_inv_powers_con_r4_vmsl, t->w_inv_powers_r4, 2 * t->n, t->q,
             VMSL_WORD_SIZE);
#endif

#ifdef AVX512_IFMA_SUPPORT
  // For avx512-ifma
  // In fact, we only need to allocate 1.25n but we allocate 2n just in case.
  _ALLOCATE_U64_ARRAY(t->w_powers_hexl, 2 * n);
  expand_w_hexl(t->w_powers_hexl, t->w_powers, n);

  _ALLOCATE_U64_ARRAY(t->w_powers_con_hexl, n * 2);
  calc_w_con(t->w_powers_con_hexl, t->w_powers_hexl, n * 2, q,
             AVX512_IFMA_WORD_SIZE);

  _ALLOCATE_U64_ARRAY(t->w_powers_r4_avx512_ifma, n * 5);
  expand_w_r4_avx512_ifma(t->w_powers_r4_avx512_ifma, t->w_powers, n, q, 0);

  _ALLOCATE_U64_ARRAY(t->w_powers_con_r4_avx512_ifma, n * 5);
  calc_w_con(t->w_powers_con_r4_avx512_ifma, t->w_powers_r4_avx512_ifma, 5 * n, q,
             AVX512_IFMA_WORD_SIZE);

  _ALLOCATE_U64_ARRAY(t->w_powers_r4_avx512_ifma_unordered, n * 5);
  expand_w_r4_avx512_ifma(t->w_powers_r4_avx512_ifma_unordered, t->w_powers, n, q,
                          1);

  _ALLOCATE_U64_ARRAY(t->w_powers_con_r4_avx512_ifma_unordered, n * 5);
  calc_w_con(t->w_powers_con_r4_avx512_ifma_unordered,
             t->w_powers_r4_avx512_ifma_unordered, n * 5, q,
             AVX512_IFMA_WORD_SIZE);

  _ALLOCATE_U64_ARRAY(t->w_powers_r4r2_avx512_ifma, n * 5);
  expand_w_r4r2_avx512_ifma(t->w_powers_r4r2_avx512_ifma, t->w_powers, n, q);

  _ALLOCATE_U64_ARRAY(t->w_powers_con_r4r2_avx512_ifma, n * 5);
  calc_w_con(t->w_powers_con_r4r2_avx512_ifma, t->w_powers_r4r2_avx512_ifma,
             n * 5, q, AVX512_IFMA_WORD_SIZE);

  _ALLOCATE_U64_ARRAY(t->w_powers_r2_16_avx512_ifma, n * 3);
  expand_w_r2_16_avx512_ifma(t->w_powers_r2_16_avx512_ifma, t->w_powers, n);

  _ALLOCATE_U64_ARRAY(t->w_powers_con_r2_16_avx512_ifma, n * 3);
  calc_w_con(t->w_powers_con_r2_16_avx512_ifma, t->w_powers_r2_16_avx512_ifma,
             n * 3, q, AVX512_IFMA_WORD_SIZE);
#endif
  return 1;
}

static inline int init_test_cases(void)
{
  for(size_t i = 0; i < NUM_OF_TEST_CASES; i++) {
    if(!_init_test(&tests[i])) {
      return 0;
    }
  }
  return 1;
}

#define FREE(p) \
  do {          \
    free((p));  \
    (p) = NULL; \
  } while(0)

static inline void _destroy_test(test_case_t *t)
{
  // for radix-2
  FREE(t->w_powers);
  FREE(t->w_powers_con);
  FREE(t->w_inv_powers);
  FREE(t->w_inv_powers_con);

  // for radix-4
  FREE(t->w_powers_r4);
  FREE(t->w_powers_con_r4);
  FREE(t->w_inv_powers_r4);
  FREE(t->w_inv_powers_con_r4);

#ifdef S390X
  // for VMSL
  FREE(t->w_powers_con_r4_vmsl);
  FREE(t->w_inv_powers_con_r4_vmsl);

#endif
#ifdef AVX512_IFMA_SUPPORT
  // for AVX512-IFMA
  FREE(t->w_powers_hexl);
  FREE(t->w_powers_con_hexl);

  FREE(t->w_powers_r4_avx512_ifma);
  FREE(t->w_powers_con_r4_avx512_ifma);

  FREE(t->w_powers_r4_avx512_ifma_unordered);
  FREE(t->w_powers_con_r4_avx512_ifma_unordered);

  FREE(t->w_powers_r4r2_avx512_ifma);
  FREE(t->w_powers_con_r4r2_avx512_ifma);

  FREE(t->w_powers_r2_16_avx512_ifma);
  FREE(t->w_powers_con_r2_16_avx512_ifma);
#endif
}

static inline void destroy_test_cases(void)
{
  for(size_t i = 0; i < NUM_OF_TEST_CASES; i++) {
    _destroy_test(&tests[i]);
  }
}

EXTERNC_END
