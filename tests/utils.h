// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

EXTERNC_BEGIN

static inline void random_buf(uint64_t *values, const size_t n, const uint64_t q)
{
  for(size_t i = 0; i < n; i++) {
    values[i] = rand() % q;
  }
}

static inline void print(UNUSED const uint64_t *values, UNUSED const size_t n)
{
#ifdef DEBUG
  for(size_t i = 0; i < n; i++) {
    printf("%lx ", values[i]);
  }
  printf("\n");
#endif
}

EXTERNC_END
