// Copyright IBM Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <stdint.h>
#include <immintrin.h>

int main(void)
{
  __m512i reg = {0};
  uint64_t mem[8] = {0};
  reg = _mm512_loadu_si512((const __m512i*)mem);
  reg = _mm512_madd52lo_epu64(reg, reg, reg);
  _mm512_storeu_si512((__m512i*)mem, reg);
  
  return 0;
}
