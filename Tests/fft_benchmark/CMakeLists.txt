cmake_minimum_required(VERSION 3.15)
project(FFT_Performance_Validation C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math")
include(FetchContent)

# Auto-detect architecture to inject the correct SIMD vectorization boundaries
if(CMAKE_SYSTEM_PROCESSOR MATCHES "amd64.*|x86_64.*|AMD64.*")
    message(STATUS "Architecture detected: x86_64. Injecting AVX2/SSE directives.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2 -mavx2")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2 -mavx2")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64.*|aarch64.*|ARM64.*")
    message(STATUS "Architecture detected: ARM64. Injecting native NEON directives.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
else()
    message(WARNING "Unknown architecture. Performance may be heavily degraded.")
endif()

find_library(FFTW3F_LIB NAMES fftw3f REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h REQUIRED)

# fftw3
add_executable(bench_fftw3 benchmark_fftw3.cpp)
target_include_directories(bench_fftw3 PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(bench_fftw3 PRIVATE ${FFTW3F_LIB} m)

# pffft
FetchContent_Declare(pffft GIT_REPOSITORY https://github.com/marton78/pffft.git GIT_TAG master)
FetchContent_MakeAvailable(pffft)
add_executable(bench_pffft benchmark_pffft.cpp)
target_link_libraries(bench_pffft PRIVATE PFFFT::PFFFT)
