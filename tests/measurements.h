// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "defs.h"

EXTERNC_BEGIN

#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#ifdef TEST_SPEED

#  ifdef INTEL_SDE
static inline void SDE_SSC_MARK(unsigned int mark_id)
{
  int ret_val;
  __asm__ __volatile__(".byte 0x64, 0x67, 0x90\n\t" // SSC mark (with ID in ebx)
                       : "=r"(ret_val)
                       : "b"(mark_id));
}

#    define SDE_SSC_START SDE_SSC_MARK(1)
#    define SDE_SSC_STOP  SDE_SSC_MARK(2)
#    define MEASURE(x) \
      SDE_SSC_START;   \
      do {             \
        x;             \
      } while(0);      \
      SDE_SSC_STOP

#  else
#    define WARMUP        10
#    define OUTER_REPEAT  10
#    define MEASURE_TIMES 200

static double start_clk;
static double end_clk;
static double total_clk;
static double temp_clk;

#    define NANO_SEC (1000000000UL)

static inline uint64_t cpucycles(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t)ts.tv_sec * NANO_SEC + ts.tv_nsec;
  ;
}

#    define MEASURE(x)                                                   \
      for(size_t warmup_itr = 0; warmup_itr < WARMUP; warmup_itr++) {    \
        {                                                                \
          x;                                                             \
        }                                                                \
      }                                                                  \
      total_clk = DBL_MAX;                                               \
      for(size_t outer_itr = 0; outer_itr < OUTER_REPEAT; outer_itr++) { \
        start_clk = cpucycles();                                         \
        for(size_t clk_itr = 0; clk_itr < MEASURE_TIMES; clk_itr++) {    \
          {                                                              \
            x;                                                           \
          }                                                              \
        }                                                                \
        end_clk  = cpucycles();                                          \
        temp_clk = (double)(end_clk - start_clk) / MEASURE_TIMES;        \
        if(total_clk > temp_clk) total_clk = temp_clk;                   \
      }                                                                  \
      printf("%9.0lu ", (uint64_t)total_clk);

#  endif
#else
#  define MEASURE(x) \
    do {             \
      x;             \
    } while(0)
#endif

EXTERNC_END
