# Copyright IBM Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "^(s390x)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DS390X")
  set(S390X 1)
elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "^(x86_64|amd64|AMD64)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DX86_64")
  set(X86_64 1)
elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "^(aarch64|arm64|arm64e)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DAARCH64")
  set(AARCH64 1)
endif()

if(X86_64)
    # Test AVX512-IFMA
    try_run(RUN_RESULT COMPILE_RESULT
            "${CMAKE_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/cmake/test_x86_64_avx512_ifma.c"
            COMPILE_DEFINITIONS "-march=native -Werror -Wall -Wpedantic"
            OUTPUT_VARIABLE OUTPUT
    )

    if(${COMPILE_RESULT} AND (RUN_RESULT EQUAL 0))
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DAVX512_IFMA_SUPPORT")
        set(AVX512_IFMA 1)
    else()
        message(STATUS "The AVX512_IFMA implementation is not supported")
    endif()
endif()
