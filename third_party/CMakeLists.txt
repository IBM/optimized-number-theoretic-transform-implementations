# Copyright IBM Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

set(THIRD_PARTY_SOURCES
    ${THIRD_PARTY_DIR}/seal/ntt_seal.c
)

if(X86_64 AND AVX512_IFMA)
    set(NTT_SOURCES ${NTT_SOURCES}
        ${SRC_DIR}/ntt_radix4_avx512_ifma.c
    )

    set(THIRD_PARTY_SOURCES ${THIRD_PARTY_SOURCES}
        ${THIRD_PARTY_DIR}/hexl/fwd-ntt-avx512.c
    )
    include_directories(${THIRD_PARTY_DIR}/hexl/)
endif()

set(TEMP ${CMAKE_C_CLANG_TIDY})
set(CMAKE_C_CLANG_TIDY "")

add_library(
  third_party OBJECT
  
  ${THIRD_PARTY_SOURCES}
)

set(CMAKE_C_CLANG_TIDY ${TEMP})
