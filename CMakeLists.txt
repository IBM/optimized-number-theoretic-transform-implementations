# Copyright IBM Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.0.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
        STRINGS "Release" "Debug" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build type (CMAKE_BUILD_TYPE): ${CMAKE_BUILD_TYPE}")

project (ntt-variants C)

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)
set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third_party)

include_directories(${INCLUDE_DIR})
include_directories(${INCLUDE_DIR}/internal)

include(cmake/arch.cmake)

include(cmake/compilation-flags.cmake)

# Depends on SRC_DIR
# and on arch.cmake
include(cmake/sources.cmake)

include(cmake/clang-format.cmake)

add_subdirectory(third_party)

add_executable(${PROJECT_NAME}
 
               ${NTT_SOURCES}
               ${MAIN_SOURCE}
               $<TARGET_OBJECTS:third_party>
)

set(BENCH ${PROJECT_NAME}-bench)
ADD_EXECUTABLE(${BENCH}

               ${NTT_SOURCES}
               ${MAIN_SOURCE}
               $<TARGET_OBJECTS:third_party>
)
SET_TARGET_PROPERTIES(${BENCH} PROPERTIES COMPILE_FLAGS "-DTEST_SPEED")
