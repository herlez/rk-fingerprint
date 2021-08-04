cmake_minimum_required(VERSION 3.15)
project(RollingHash)

set(CMAKE_CXX_STANDARD 20)
add_executable(bench benchmark.cpp)
target_compile_options(bench PRIVATE -Wall -Wextra -pedantic -Werror -O3 -save-temps)
