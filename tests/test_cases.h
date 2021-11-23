// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdlib.h>

#include "fast_mul_operators.h"
#include "pre_compute.h"

EXTERNC_BEGIN

/*
 * The pointers in this file are 64 bytes aligned.
 */
typedef struct aligned64_ptr_s {
  void *    base;
  uint64_t *ptr;
} aligned64_ptr_t;

typedef struct unaligned64_ptr_s {
  void *    base;
  uint64_t *ptr;
} unaligned64_ptr_t;

static inline int allocate_aligned_array(aligned64_ptr_t *aptr, size_t qw_num)
{
  size_t size_to_allocate = qw_num * sizeof(uint64_t) + 64;
  if(NULL == ((aptr->base) = malloc(size_to_allocate))) {
    printf("Allocation error");
    return ERROR;
  }
  aptr->ptr = (uint64_t *)(((uint64_t)aptr->base & (~0x3fULL)) + 64);
  return SUCCESS;
}

static inline int allocate_unaligned_array(unaligned64_ptr_t *aptr, size_t qw_num)
{
  size_t size_to_allocate = qw_num * sizeof(uint64_t) + 64 + 8;
  if(NULL == ((aptr->base) = malloc(size_to_allocate))) {
    printf("Allocation error");
    return ERROR;
  }
  aptr->ptr = (uint64_t *)(((uint64_t)aptr->base & (~0x3fULL)) + 64 + 8);
  return SUCCESS;
}

static inline void free_aligned_array(aligned64_ptr_t *aptr)
{
  free(aptr->base);
  aptr->base = NULL;
  aptr->ptr  = NULL;
}

static inline void free_unaligned_array(unaligned64_ptr_t *aptr)
{
  free(aptr->base);
  aptr->base = NULL;
  aptr->ptr  = NULL;
}

typedef struct test_case_s {
  // These parameters are predefined
  uint64_t m;
  uint64_t q;
  uint64_t w;
  uint64_t w_inv; // w^(-1) mod q
  mul_op_t n_inv; // 2^(-m) mod q

  // These parameters are dinamically computed based on the above values.
  uint64_t        n;
  uint64_t        qneg;
  uint64_t        q2;
  uint64_t        q4;
  aligned64_ptr_t w_powers;
  aligned64_ptr_t w_powers_con;
  aligned64_ptr_t w_inv_powers;
  aligned64_ptr_t w_inv_powers_con;

  // For radix-4 tests
  aligned64_ptr_t w_powers_r4;
  aligned64_ptr_t w_powers_con_r4;
  aligned64_ptr_t w_inv_powers_r4;
  aligned64_ptr_t w_inv_powers_con_r4;

#ifdef S390X
  // For radix-4 tests with VMSL (56-bits instead of 64-bits)
  aligned64_ptr_t w_powers_con_r4_vmsl;
  aligned64_ptr_t w_inv_powers_con_r4_vmsl;
  mul_op_t        n_inv_vmsl;
#endif

#ifdef AVX512_IFMA_SUPPORT
  // For radix-2 tests with AVX512-IFMA on X86-64 bit platofrms (52-bits)
  aligned64_ptr_t w_powers_hexl;
  aligned64_ptr_t w_powers_con_hexl;

  aligned64_ptr_t w_powers_r4_avx512_ifma;
  aligned64_ptr_t w_powers_con_r4_avx512_ifma;

  aligned64_ptr_t w_powers_r4_avx512_ifma_unordered;
  aligned64_ptr_t w_powers_con_r4_avx512_ifma_unordered;

  aligned64_ptr_t w_powers_r4r2_avx512_ifma;
  aligned64_ptr_t w_powers_con_r4r2_avx512_ifma;

  aligned64_ptr_t w_powers_r2_16_avx512_ifma;
  aligned64_ptr_t w_powers_con_r2_16_avx512_ifma;
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
   .q        = 0x10001, // NOLINT
   .w        = 3,
   .w_inv    = 21846,  // NOLINT
   .n_inv.op = 65535}, // NOLINT
  {.m        = 15,
   .q        = 0x80000001c0001, // NOLINT
   .w        = 82138512871,
   .w_inv    = 535648572761016,   // NOLINT
   .n_inv.op = 2251731096043465}, // NOLINT
  {.m        = 16,                // NOLINT
   .q        = 0x7ffe0001,        // NOLINT
   .w        = 1859,
   .w_inv    = 1579037640, // NOLINT
   .n_inv.op = 2147319811},
  {.m        = 16,              // NOLINT
   .q        = 0x7fffffffe0001, // NOLINT
   .w        = 29454831443,
   .w_inv    = 520731633805630,   // NOLINT
   .n_inv.op = 2251765453815811}, // NOLINT
  {.m        = 17,                // NOLINT
   .q        = 0x100180001,       // NOLINT
   .w        = 79247,
   .w_inv    = 4203069932,   // NOLINT
   .n_inv.op = 4296507381}}; // NOLINT

#define NUM_OF_TEST_CASES (sizeof(tests) / sizeof(test_case_t))

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
  allocate_aligned_array(&t->w_powers, n);
  calc_w(t->w_powers.ptr, w, n, q, m);

  allocate_aligned_array(&t->w_powers_con, n);
  calc_w_con(t->w_powers_con.ptr, t->w_powers.ptr, n, q, WORD_SIZE);

  allocate_aligned_array(&t->w_inv_powers, n);
  calc_w_inv(t->w_inv_powers.ptr, w_inv, n, q, m);

  allocate_aligned_array(&t->w_inv_powers_con, n);
  calc_w_con(t->w_inv_powers_con.ptr, t->w_inv_powers.ptr, n, q, WORD_SIZE);

  // Expand the list of powers to support the radix-4 case.
  allocate_aligned_array(&t->w_powers_r4, 2 * n);
  expand_w(t->w_powers_r4.ptr, t->w_powers.ptr, n, q);

  allocate_aligned_array(&t->w_powers_con_r4, 2 * n);
  calc_w_con(t->w_powers_con_r4.ptr, t->w_powers_r4.ptr, 2 * n, q, WORD_SIZE);

  allocate_aligned_array(&t->w_inv_powers_r4, 2 * n);
  expand_w(t->w_inv_powers_r4.ptr, t->w_inv_powers.ptr, n, q);

  allocate_aligned_array(&t->w_inv_powers_con_r4, 2 * n);
  calc_w_con(t->w_inv_powers_con_r4.ptr, t->w_inv_powers_r4.ptr, 2 * n, q,
             WORD_SIZE);

#ifdef S390X
  t->n_inv_vmsl.con = calc_ninv_con(&t->n_inv.op, q, VMSL_WORD_SIZE);
  t->n_inv_vmsl.op  = t->n_inv.op;

  // for radix-4 vmsl
  allocate_aligned_array(&t->w_powers_con_r4_vmsl, 2 * n);
  calc_w_con(t->w_powers_con_r4_vmsl.ptr, t->w_powers_r4.ptr, 2 * t->n, t->q,
             VMSL_WORD_SIZE);

  allocate_aligned_array(&t->w_inv_powers_con_r4_vmsl, 2 * n);
  calc_w_con(t->w_inv_powers_con_r4_vmsl.ptr, t->w_inv_powers_r4.ptr, 2 * t->n,
             t->q, VMSL_WORD_SIZE);
#endif

#ifdef AVX512_IFMA_SUPPORT
  // For avx512-ifma
  // In fact, we only need to allocate 1.25n but we allocate 2n just in case.
  allocate_aligned_array(&t->w_powers_hexl, 2 * n);
  expand_w_hexl(t->w_powers_hexl.ptr, t->w_powers.ptr, n);

  allocate_aligned_array(&t->w_powers_con_hexl, n * 2);
  calc_w_con(t->w_powers_con_hexl.ptr, t->w_powers_hexl.ptr, n * 2, q,
             AVX512_IFMA_WORD_SIZE);

  allocate_aligned_array(&t->w_powers_r4_avx512_ifma, n * 5);
  expand_w_r4_avx512_ifma(t->w_powers_r4_avx512_ifma.ptr, t->w_powers.ptr, n, q,
                          0);

  allocate_aligned_array(&t->w_powers_con_r4_avx512_ifma, n * 5);
  calc_w_con(t->w_powers_con_r4_avx512_ifma.ptr, t->w_powers_r4_avx512_ifma.ptr,
             5 * n, q, AVX512_IFMA_WORD_SIZE);

  allocate_aligned_array(&t->w_powers_r4_avx512_ifma_unordered, n * 5);
  expand_w_r4_avx512_ifma(t->w_powers_r4_avx512_ifma_unordered.ptr,
                          t->w_powers.ptr, n, q, 1);

  allocate_aligned_array(&t->w_powers_con_r4_avx512_ifma_unordered, n * 5);
  calc_w_con(t->w_powers_con_r4_avx512_ifma_unordered.ptr,
             t->w_powers_r4_avx512_ifma_unordered.ptr, n * 5, q,
             AVX512_IFMA_WORD_SIZE);

  allocate_aligned_array(&t->w_powers_r4r2_avx512_ifma, n * 5);
  expand_w_r4r2_avx512_ifma(t->w_powers_r4r2_avx512_ifma.ptr, t->w_powers.ptr, n,
                            q);

  allocate_aligned_array(&t->w_powers_con_r4r2_avx512_ifma, n * 5);
  calc_w_con(t->w_powers_con_r4r2_avx512_ifma.ptr,
             t->w_powers_r4r2_avx512_ifma.ptr, n * 5, q, AVX512_IFMA_WORD_SIZE);

  allocate_aligned_array(&t->w_powers_r2_16_avx512_ifma, n * 3);
  expand_w_r2_16_avx512_ifma(t->w_powers_r2_16_avx512_ifma.ptr, t->w_powers.ptr,
                             n);

  allocate_aligned_array(&t->w_powers_con_r2_16_avx512_ifma, n * 3);
  calc_w_con(t->w_powers_con_r2_16_avx512_ifma.ptr,
             t->w_powers_r2_16_avx512_ifma.ptr, n * 3, q, AVX512_IFMA_WORD_SIZE);
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

static inline void _destroy_test(test_case_t *t)
{
  // for radix-2
  free_aligned_array(&t->w_powers);
  free_aligned_array(&t->w_powers_con);
  free_aligned_array(&t->w_inv_powers);
  free_aligned_array(&t->w_inv_powers_con);

  // for radix-4
  free_aligned_array(&t->w_powers_r4);
  free_aligned_array(&t->w_powers_con_r4);
  free_aligned_array(&t->w_inv_powers_r4);
  free_aligned_array(&t->w_inv_powers_con_r4);

#ifdef S390X
  // for VMSL
  free_aligned_array(&t->w_powers_con_r4_vmsl);
  free_aligned_array(&t->w_inv_powers_con_r4_vmsl);

#endif
#ifdef AVX512_IFMA_SUPPORT
  // for AVX512-IFMA
  free_aligned_array(&t->w_powers_hexl);
  free_aligned_array(&t->w_powers_con_hexl);

  free_aligned_array(&t->w_powers_r4_avx512_ifma);
  free_aligned_array(&t->w_powers_con_r4_avx512_ifma);

  free_aligned_array(&t->w_powers_r4_avx512_ifma_unordered);
  free_aligned_array(&t->w_powers_con_r4_avx512_ifma_unordered);

  free_aligned_array(&t->w_powers_r4r2_avx512_ifma);
  free_aligned_array(&t->w_powers_con_r4r2_avx512_ifma);

  free_aligned_array(&t->w_powers_r2_16_avx512_ifma);
  free_aligned_array(&t->w_powers_con_r2_16_avx512_ifma);
#endif
}

static inline void destroy_test_cases(void)
{
  for(size_t i = 0; i < NUM_OF_TEST_CASES; i++) {
    _destroy_test(&tests[i]);
  }
}

EXTERNC_END
