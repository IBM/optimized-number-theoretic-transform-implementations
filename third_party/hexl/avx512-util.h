/*
 * The code in this file was extracted and modified from 
b) Intel HEXL GitHub (fwd-ntt-avx512.cpp)[https://github.com/intel/hexl/blob/db9535c140227010c5c9d6f34a11054b16f02de7/hexl/ntt/fwd-ntt-avx512.cpp] commit 4d9806f 01 Sep 2021 that has an (Apache 2.0 license)[https://github.com/intel/hexl/blob/main/LICENSE]

The code was converted from C++ to native C by
- Removing templates and converting all relevant functions to `static inline` functions.
- Converting `reinterpret_cast` and `static_cast` to C-style casting.
- Fixed the `BitShift` parameter to 52 and removed code paths (and branches) to its other values.
- Set the `InputLessThanMod` parameter as an input parameter to the relevant functions.
- Converted the `HEXL_LOOP_UNROLL_N` macros to `LOOP_UNROLL_N` macros.
- Defined the `HEXL_CHECK` and ``HEXL_VLOG` macros as empty macros.
*/

#pragma once

#include <immintrin.h>

#define HEXL_CHECK(...) {}
#define HEXL_CHECK_BOUNDS(...) {}
#define HEXL_VLOG(...) {}

static inline __m512i _mm512_hexl_mullo_epi_64(__m512i x, __m512i y) {
  return _mm512_mullo_epi64(x, y);
}

static inline __m512i _mm512_hexl_mullo_epi_52(__m512i x, __m512i y) {
  __m512i zero = _mm512_set1_epi64(0);
  return _mm512_madd52lo_epu64(zero, x, y);
}

static inline __m512i _mm512_hexl_mullo_add_lo_epi_52(__m512i x, __m512i y,
                                                __m512i z) {
  __m512i result = _mm512_madd52lo_epu64(x, y, z);

  // Clear high 12 bits from result
  const __m512i two_pow52_min1 = _mm512_set1_epi64((1ULL << 52) - 1);
  result = _mm512_and_epi64(result, two_pow52_min1);
  return result;
}

static inline __m512i _mm512_hexl_mullo_add_lo_epi_64(__m512i x, __m512i y,
                                                __m512i z) {
  __m512i prod = _mm512_mullo_epi64(y, z);
  return _mm512_add_epi64(x, prod);
}

// Returns x mod q across each 64-bit integer SIMD lanes
// Assumes x < InputModFactor * q in all lanes
static const int InputModFactor = 2;
static inline __m512i _mm512_hexl_small_mod_epu64(__m512i x, __m512i q) {
  __m512i* q_times_2 = NULL;
  __m512i* q_times_4 = NULL;
  HEXL_CHECK(InputModFactor == 1 || InputModFactor == 2 ||
                 InputModFactor == 4 || InputModFactor == 8,
             "InputModFactor must be 1, 2, 4, or 8");
  if (InputModFactor == 1) {
    return x;
  }
  if (InputModFactor == 2) {
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
  }
  if (InputModFactor == 4) {
    HEXL_CHECK(q_times_2 != nullptr, "q_times_2 must not be nullptr");
    x = _mm512_min_epu64(x, _mm512_sub_epi64(x, *q_times_2));
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
  }
  if (InputModFactor == 8) {
    HEXL_CHECK(q_times_2 != nullptr, "q_times_2 must not be nullptr");
    HEXL_CHECK(q_times_4 != nullptr, "q_times_4 must not be nullptr");
    x = _mm512_min_epu64(x, _mm512_sub_epi64(x, *q_times_4));
    x = _mm512_min_epu64(x, _mm512_sub_epi64(x, *q_times_2));
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
  }
  HEXL_CHECK(false, "Invalid InputModFactor");
  return x;  // Return dummy value
}

static inline __m512i _mm512_hexl_mulhi_epi_64(__m512i x, __m512i y) {
  // https://stackoverflow.com/questions/28807341/simd-signed-with-unsigned-multiplication-for-64-bit-64-bit-to-128-bit
  __m512i lo_mask = _mm512_set1_epi64(0x00000000ffffffff);
  // Shuffle high bits with low bits in each 64-bit integer =>
  // x0_lo, x0_hi, x1_lo, x1_hi, x2_lo, x2_hi, ...
  __m512i x_hi = _mm512_shuffle_epi32(x, (_MM_PERM_ENUM)0xB1);
  // y0_lo, y0_hi, y1_lo, y1_hi, y2_lo, y2_hi, ...
  __m512i y_hi = _mm512_shuffle_epi32(y, (_MM_PERM_ENUM)0xB1);
  __m512i z_lo_lo = _mm512_mul_epu32(x, y);        // x_lo * y_lo
  __m512i z_lo_hi = _mm512_mul_epu32(x, y_hi);     // x_lo * y_hi
  __m512i z_hi_lo = _mm512_mul_epu32(x_hi, y);     // x_hi * y_lo
  __m512i z_hi_hi = _mm512_mul_epu32(x_hi, y_hi);  // x_hi * y_hi

  //                   x_hi | x_lo
  // x                 y_hi | y_lo
  // ------------------------------
  //                  [x_lo * y_lo]    // z_lo_lo
  // +           [z_lo * y_hi]         // z_lo_hi
  // +           [x_hi * y_lo]         // z_hi_lo
  // +    [x_hi * y_hi]                // z_hi_hi
  //     ^-----------^ <-- only bits needed
  //  sum_|  hi | mid | lo  |

  // Low bits of z_lo_lo are not needed
  __m512i z_lo_lo_shift = _mm512_srli_epi64(z_lo_lo, 32);

  //                   [x_lo  *  y_lo] // z_lo_lo
  //          + [z_lo  *  y_hi]        // z_lo_hi
  //          ------------------------
  //            |    sum_tmp   |
  //            |sum_mid|sum_lo|
  __m512i sum_tmp = _mm512_add_epi64(z_lo_hi, z_lo_lo_shift);
  __m512i sum_lo = _mm512_and_si512(sum_tmp, lo_mask);
  __m512i sum_mid = _mm512_srli_epi64(sum_tmp, 32);
  //            |       |sum_lo|
  //          + [x_hi   *  y_lo]       // z_hi_lo
  //          ------------------
  //            [   sum_mid2   ]
  __m512i sum_mid2 = _mm512_add_epi64(z_hi_lo, sum_lo);
  __m512i sum_mid2_hi = _mm512_srli_epi64(sum_mid2, 32);
  __m512i sum_hi = _mm512_add_epi64(z_hi_hi, sum_mid);
  return _mm512_add_epi64(sum_hi, sum_mid2_hi);
}

static inline __m512i _mm512_hexl_mulhi_approx_epi_64(__m512i x, __m512i y) {
  // https://stackoverflow.com/questions/28807341/simd-signed-with-unsigned-multiplication-for-64-bit-64-bit-to-128-bit
  __m512i lo_mask = _mm512_set1_epi64(0x00000000ffffffff);
  // Shuffle high bits with low bits in each 64-bit integer =>
  // x0_lo, x0_hi, x1_lo, x1_hi, x2_lo, x2_hi, ...
  __m512i x_hi = _mm512_shuffle_epi32(x, (_MM_PERM_ENUM)0xB1);
  // y0_lo, y0_hi, y1_lo, y1_hi, y2_lo, y2_hi, ...
  __m512i y_hi = _mm512_shuffle_epi32(y, (_MM_PERM_ENUM)0xB1);
  __m512i z_lo_hi = _mm512_mul_epu32(x, y_hi);     // x_lo * y_hi
  __m512i z_hi_lo = _mm512_mul_epu32(x_hi, y);     // x_hi * y_lo
  __m512i z_hi_hi = _mm512_mul_epu32(x_hi, y_hi);  // x_hi * y_hi

  //                   x_hi | x_lo
  // x                 y_hi | y_lo
  // ------------------------------
  //                  [x_lo * y_lo]    // unused, resulting in approximation
  // +           [z_lo * y_hi]         // z_lo_hi
  // +           [x_hi * y_lo]         // z_hi_lo
  // +    [x_hi * y_hi]                // z_hi_hi
  //     ^-----------^ <-- only bits needed
  //  sum_|  hi | mid | lo  |

  __m512i sum_lo = _mm512_and_si512(z_lo_hi, lo_mask);
  __m512i sum_mid = _mm512_srli_epi64(z_lo_hi, 32);
  //            |       |sum_lo|
  //          + [x_hi   *  y_lo]       // z_hi_lo
  //          ------------------
  //            [   sum_mid2   ]
  __m512i sum_mid2 = _mm512_add_epi64(z_hi_lo, sum_lo);
  __m512i sum_mid2_hi = _mm512_srli_epi64(sum_mid2, 32);
  __m512i sum_hi = _mm512_add_epi64(z_hi_hi, sum_mid);
  return _mm512_add_epi64(sum_hi, sum_mid2_hi);
}

static inline __m512i _mm512_hexl_mulhi_approx_epi_52(__m512i x, __m512i y) {
  __m512i zero = _mm512_set1_epi64(0);
  return _mm512_madd52hi_epu64(zero, x, y);
}

static inline __m512i _mm512_hexl_mulhi_epi_52(__m512i x, __m512i y) {
  __m512i zero = _mm512_set1_epi64(0);
  return _mm512_madd52hi_epu64(zero, x, y);
}
