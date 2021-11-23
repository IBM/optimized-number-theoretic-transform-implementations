// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "pre_compute.h"
#include "tests.h"

int main(UNUSED int argc, UNUSED char *argv[])
{
  init_test_cases();

#ifdef TEST_SPEED
  if(argc == 2) {
    printf("Testing test 9 and func %ld cycle=", strtol(argv[1], NULL, 0));
    test_fwd_single_case(&tests[9], strtol(argv[1], NULL, 0));
    printf("\n");
    return SUCCESS;
  }

  printf("\n\nTesting forward NTT with unaligned inputs\n\n");
  report_test_fwd_perf_headers();
  for(size_t i = 0; i < NUM_OF_TEST_CASES; i++) {
    test_unaligned_fwd_perf(&tests[i]);
  }

  printf("Testing forward NTT with aligned inputs\n\n");
  report_test_fwd_perf_headers();
  for(size_t i = 0; i < NUM_OF_TEST_CASES; i++) {
    test_aligned_fwd_perf(&tests[i]);
  }

  printf("Testing inverse NTT with unaligned inputs\n\n");
  report_test_inv_perf_headers();
  for(size_t i = 0; i < NUM_OF_TEST_CASES; i++) {
    test_inv_perf(&tests[i]);
  }

#else

  for(size_t i = 0; i < NUM_OF_TEST_CASES; i++) {
    printf("Test %2.0lu\n", i);
    if(SUCCESS != test_correctness(&tests[i])) {
      destroy_test_cases();
      return SUCCESS;
    }
  }
#endif

  destroy_test_cases();
  return SUCCESS;
}
