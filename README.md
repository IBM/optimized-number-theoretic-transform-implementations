# optimized-number-theoretic-transform-implementations

This sample code package is an implementation of the Number Theoretic Transform (NTT) 
algorithm for the ring R/(X^N + 1) where N=2^m.

This sample code provides testing binaries but no shared or static libraries. This is because the code is designed to be used for benchmarking purposes only and not in final products.

## License

This project is licensed under the Apache-2.0 License.

Dependencies
-----
This package requires 
- CMake 3 and above 
- A compiler that supports the required C intrinsics (e.g., VMSL on s390x machines or AVX512-IFMA on X86_64 platforms). For example, GCC-10 and Clang-12.

BUILD
-----

To build the directory first create a working directory
```
mkdir build
cd build
```

Then, run CMake and compile
```
cmake ..
make
```

To run

`./ntt-variants`

Additional CMake compilation flags:
  - DEBUG       - To enable debug prints

To clean - remove the `build` directory. Note that a "clean" is required prior to compilation with modified flags.

To format (`clang-format-9` or above is required):

`make format`

To use clang-tidy (`clang-tidy-9` is required):

```
CC=clang-12 cmake -DCMAKE_C_CLANG_TIDY="clang-tidy;--format-style=file" ..
make 
```

Before committing code, please test it using
`tests/pre-commit-script.sh` 
This will run all the sanitizers and also `clang-format` and `clang-tidy` (requires clang-9 to be installed).

The package was compiled and tested with gcc-10 and clang-12 in 64-bit mode. 
Tests were run on a Linux (Ubuntu 20.04.2 LTS) OS on s390x and on ICL x86-64 platforms. 
Compilation on other platforms may require some adjustments.

Performance measurements
------------------------
The performance measurements are reported in processor cycles (per single core). The results are obtained using the following methodology. Each measured function was isolated, run 10 times (warm-up), followed by 200 iterations that were clocked and averaged. To minimize the effect of background tasks running on the system, every experiment was repeated 10 times, and the minimum result is reported.

To run the benchmarking

`./ntt-variants-bench`

Testing
-------
- The library has several fixed test-cases with different values of `q` and `N`. 
- The library was run using Address/Undefined-Behaviour sanitizers.

Notes
-----
clang-12 achieves much better results than GCC-10.
