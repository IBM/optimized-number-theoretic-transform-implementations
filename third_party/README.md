Third party code
----------------
The code in this directory was extracted and modified from 
1) Microsoft SEAL GitHub file (dwthandler.h)[https://github.com/microsoft/SEAL/blob/d045f1beff96dff0fccc7fa0c5acb1493a65338c/native/src/seal/util/dwthandler.h] commit d045f1b on 15 Jun 2021 that has an (MIT license)[https://github.com/microsoft/SEAL/blob/main/LICENSE]
2) Intel HEXL GitHub (fwd-ntt-avx512.cpp)[https://github.com/intel/hexl/blob/db9535c140227010c5c9d6f34a11054b16f02de7/hexl/ntt/fwd-ntt-avx512.cpp] commit 4d9806f 01 Sep 2021 that has an (Apache 2.0 license)[https://github.com/intel/hexl/blob/main/LICENSE]

The code was converted from C++ to native C by
- Removing templates and converting all relevant functions to `static inline` functions.
- Converting `reinterpret_cast` and `static_cast` to C-style casting.

Specifically for HEXL we
- Fixed the `BitShift` parameter to 52 and removed code paths (and branches) to its other values.
- Set the `InputLessThanMod` parameter as an input parameter to the relevant functions.
- Converted the `HEXL_LOOP_UNROLL_N` macros to `LOOP_UNROLL_N` macros.
- Defined the `HEXL_CHECK` and ``HEXL_VLOG` macros as empty macros.
