// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

EXTERNC_BEGIN

#include "test_cases.h"

typedef enum
{
  FIRST_FWD = 0,
  FWD_REF   = FIRST_FWD,
  FWD_SEAL,
  FWD_R4,
  FWD_R4x4,
  FWD_R4_VMSL,
  FWD_R4_HEXL,
  FWD_R4_AVX512_IFMA,
  FWD_R4_AVX512_IFMA_UNORDERED,
  FWD_R4R2_AVX512_IFMA,
  FWD_R2_R16_AVX512_IFMA,
  MAX_FWD = FWD_R2_R16_AVX512_IFMA
} func_num_t;

#ifdef TEST_SPEED

void report_test_fwd_perf_headers(void);
void report_test_inv_perf_headers(void);

void test_fwd_perf(const test_case_t *t);
void test_inv_perf(const test_case_t *t);

void test_fwd_single_case(const test_case_t *t, func_num_t func_num);

#else

int test_correctness(const test_case_t *t);

#endif

EXTERNC_END
