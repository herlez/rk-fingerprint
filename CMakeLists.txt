cmake_minimum_required(VERSION 3.15)
project(RollingHash)

set(CMAKE_CXX_STANDARD 20)
add_executable(bench benchmark.cpp)
target_compile_options(bench PRIVATE -Wall -Wextra -pedantic -O3 -save-temps) # -Winline
target_include_directories(bench PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extlib/kr-fingerprinting>
  $<INSTALL_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extlib/kr-fingerprinting>
)
