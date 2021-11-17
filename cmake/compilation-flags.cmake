# Copyright IBM Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(CLANG 1)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb -O3 -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden -Wall -Wextra -Werror -Wpedantic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wunused -Wcomment -Wchar-subscripts -Wuninitialized -Wshadow")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wwrite-strings -Wformat-security -Wcast-qual -Wunused-result")

# Remove any occurrences of parameter -O from CMAKE_C_FLAGS_<BuildType>
string(TOUPPER "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}" varname)
string(REGEX REPLACE "-O[0-3]?" "" ${varname} "${${varname}}")
unset(varname)

if(S390X)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=z14 -mvx -mzvector")
    if (NOT CLANG)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mbranch-cost=3")
    endif()
elseif(X86_64)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mno-red-zone")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native")
endif()

if(MSAN)
    if(NOT CLANG)
        message(FATAL_ERROR "Cannot enable MSAN unless using Clang")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer")
endif()

if(ASAN)
    if(NOT CLANG)
        message(FATAL_ERROR "Cannot enable ASAN unless using Clang")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fsanitize-address-use-after-scope -fno-omit-frame-pointer")
endif()

if(TSAN)
    if(NOT CLANG)
        message(FATAL_ERROR "Cannot enable TSAN unless using Clang")
    endif()
    if(S390X)
        message(FATAL_ERROR "Cannot enable TSAN for s390x machines")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread")
endif()

if(UBSAN)
    if(NOT CLANG)
        message(FATAL_ERROR "Cannot enable UBSAN unless using Clang")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined")
endif()

if(DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEBUG")
endif()

if(INTEL_SDE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINTEL_SDE")
endif()
