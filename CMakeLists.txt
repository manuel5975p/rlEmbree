cmake_minimum_required(VERSION 3.20)
project(rlembree)
include(FetchContent)
include(CheckCXXSourceCompiles)

# Check for AVX support
set(CMAKE_REQUIRED_FLAGS "-march=native")
check_cxx_source_compiles("
#include <immintrin.h>
int main() {
    __m128i a;
    _mm_test_all_ones(a);
    return 0;
}" HAS_SSE42_SUPPORT)

check_cxx_source_compiles("
#include <immintrin.h>
int main() {
    __m256 a = _mm256_setzero_ps();
    return 0;
}" HAS_AVX_SUPPORT)

# Check for AVX2 support
check_cxx_source_compiles("
#include <immintrin.h>
int main() {
    __m256i a = _mm256_setzero_si256();
    return 0;
}" HAS_AVX2_SUPPORT)

# Check for AVX512 support
check_cxx_source_compiles("
#include <immintrin.h>
int main() {
    __m512i a = _mm512_setzero_si512();
    return 0;
}" HAS_AVX512_SUPPORT)
set(MAXSUP "SSE2")
# Report results
if(HAS_SSE42_SUPPORT)
    message(STATUS "SSE42 is supported.")
    set(MAXSUP "SSE42")
else()
    message(STATUS "SSE42 is not supported.")
endif()
if(HAS_AVX_SUPPORT)
    message(STATUS "AVX is supported.")
    set(MAXSUP "AVX")
else()
    message(STATUS "AVX is not supported.")
endif()

if(HAS_AVX2_SUPPORT)
    message(STATUS "AVX2 is supported.")
    set(MAXSUP "AVX2")
else()
    message(STATUS "AVX2 is not supported.")
endif()

if(HAS_AVX512_SUPPORT)
    message(STATUS "AVX512 is supported.")
    set(MAXSUP "AVX512")
else()
    message(STATUS "AVX512 is not supported.")
endif()
set(EMBREE_GEOMETRY_CURVE        OFF CACHE BOOL "" FORCE)
set(EMBREE_GEOMETRY_GRID         OFF CACHE BOOL "" FORCE)
set(EMBREE_GEOMETRY_USER         OFF CACHE BOOL "" FORCE)
set(EMBREE_GEOMETRY_SUBDIVISION  OFF CACHE BOOL "" FORCE)
set(EMBREE_ISA_AVX512            OFF CACHE BOOL "" FORCE)
set(EMBREE_ISA_AVX2              OFF CACHE BOOL "" FORCE)
set(EMBREE_ISA_AVX               OFF CACHE BOOL "" FORCE)
set(EMBREE_ISA_SSE2              OFF CACHE BOOL "" FORCE)
set(EMBREE_ISA_SSE42             OFF CACHE BOOL "" FORCE)
set(EMBREE_GEOMETRY_USER         OFF CACHE BOOL "" FORCE)
set(EMBREE_GEOMETRY_USER         OFF CACHE BOOL "" FORCE)
set(EMBREE_GEOMETRY_USER         OFF CACHE BOOL "" FORCE)
set(EMBREE_TUTORIALS             OFF CACHE BOOL "" FORCE)
set(EMBREE_TASKING_SYSTEM   INTERNAL CACHE STRING "" FORCE)
set(EMBREE_TUTORIALS_GLFW        OFF CACHE BOOL "" FORCE)
#set(EMBREE_STATIC_LIB            ON  CACHE BOOL "" FORCE)
#set(EMBREE_RAY_PACKETS           OFF CACHE BOOL "" FORCE)
set("EMBREE_ISA_${MAXSUP}"       ON  CACHE BOOL "" FORCE)
FetchContent_Declare(embreec 
    URL https://github.com/manuel5975p/manuel5975p.github.io/raw/master/slimbree-4.3.1.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable(embreec)

add_executable(conversion conversion.cpp)
target_compile_features(conversion PUBLIC cxx_std_20)
target_compile_options(conversion PUBLIC "-fopenmp")
target_link_libraries(conversion embree raylib)
target_link_options(conversion PUBLIC "-fopenmp")

add_executable         (cfile cfile.c conversion.cpp)
target_compile_features(cfile PUBLIC c_std_11)
target_compile_options (cfile PUBLIC "-fopenmp")
target_link_libraries  (cfile embree raylib)
target_link_options    (cfile PUBLIC "-fopenmp")