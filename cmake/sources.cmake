# Copyright IBM Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

set(NTT_SOURCES 
    ${SRC_DIR}/ntt_radix4.c
    ${SRC_DIR}/ntt_radix4x4.c
    ${SRC_DIR}/ntt_reference.c
)

if(S390X)
    set(NTT_SOURCES ${NTT_SOURCES}
        ${SRC_DIR}/ntt_radix4_s390x_vef.c
    )
endif()

if(X86_64 AND AVX512_IFMA)
    set(NTT_SOURCES ${NTT_SOURCES}
        ${SRC_DIR}/ntt_radix4_avx512_ifma.c
        ${SRC_DIR}/ntt_r4r2_avx512_ifma.c
        ${SRC_DIR}/ntt_r2_16_avx512_ifma.c
        ${SRC_DIR}/ntt_radix4_avx512_ifma_unordered.c
    )
endif()

set(MAIN_SOURCE 
    ${TESTS_DIR}/main.c
    ${TESTS_DIR}/bench.c
    ${TESTS_DIR}/test_correctness.c
)
