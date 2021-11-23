// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
#  define EXTERNC       extern "C"
#  define EXTERNC_BEGIN extern "C" {
#  define EXTERNC_END   }
#else
#  define EXTERNC
#  define EXTERNC_BEGIN
#  define EXTERNC_END
#endif

#define SUCCESS 0
#define ERROR   (-1)

#define GUARD(func)         \
  {                         \
    if(SUCCESS != (func)) { \
      return ERROR;         \
    }                       \
  }

#define GUARD_MSG(func, msg) \
  {                          \
    if(SUCCESS != (func)) {  \
      printf(msg);           \
      return ERROR;          \
    }                        \
  }

#if defined(__GNUC__) || defined(__clang__)
#  define UNUSED __attribute__((unused))
#else
#  define UNUSED
#endif

#define WORD_SIZE             64UL
#define VMSL_WORD_SIZE        56UL
#define AVX512_IFMA_WORD_SIZE 52UL

#if WORD_SIZE == 64
#  define WORD_SIZE_MASK (-1UL)
#else
#  define WORD_SIZE_MASK ((1UL << WORD_SIZE) - 1)
#endif

#define HIGH_WORD(x) ((x) >> WORD_SIZE)
#define LOW_WORD(x)  ((x)&WORD_SIZE_MASK)

#define VMSL_WORD_SIZE_MASK ((1UL << VMSL_WORD_SIZE) - 1)
#define HIGH_VMSL_WORD(x)   (uint64_t)((__uint128_t)(x) >> VMSL_WORD_SIZE)
#define LOW_VMSL_WORD(x)    ((x)&VMSL_WORD_SIZE_MASK)

#define AVX512_IFMA_WORD_SIZE_MASK   ((1UL << AVX512_IFMA_WORD_SIZE) - 1)
#define AVX512_IFMA_MAX_MODULUS      49UL
#define AVX512_IFMA_MAX_MODULUS_MASK (~((1UL << AVX512_IFMA_MAX_MODULUS) - 1))

// Check whether N=2^m where m is odd by masking it.
#define ODD_POWER_MASK  0xaaaaaaaaaaaaaaaa
#define REM1_POWER_MASK 0x2222222222222222
#define REM2_POWER_MASK 0x4444444444444444
#define REM3_POWER_MASK 0x8888888888888888

#define HAS_AN_EVEN_POWER(n) (!((n)&ODD_POWER_MASK))
#define HAS_AN_REM1_POWER(n) ((n)&REM1_POWER_MASK)
#define HAS_AN_REM2_POWER(n) ((n)&REM2_POWER_MASK)
#define HAS_AN_REM3_POWER(n) ((n)&REM3_POWER_MASK)

#if defined(__GNUC__) && (__GNUC__ >= 8)
#  define GCC_SUPPORT_UNROLL_PRAGMA
#endif

#ifdef GCC_SUPPORT_UNROLL_PRAGMA
#  define LOOP_UNROLL_2 _Pragma("GCC unroll 2")
#  define LOOP_UNROLL_4 _Pragma("GCC unroll 4")
#  define LOOP_UNROLL_8 _Pragma("GCC unroll 8")
#elif defined(__clang__)
#  define LOOP_UNROLL_2 _Pragma("clang loop unroll_count(2)")
#  define LOOP_UNROLL_4 _Pragma("clang loop unroll_count(4)")
#  define LOOP_UNROLL_8 _Pragma("clang loop unroll_count(8)")
#else
#  define LOOP_UNROLL_2
#  define LOOP_UNROLL_4
#  define LOOP_UNROLL_8
#endif

#define ALIGN(n) __attribute__((aligned(n)))
