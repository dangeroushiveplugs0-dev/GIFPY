cmake_minimum_required(VERSION 3.20)
project(GIFPY LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(gifpy_liquid src/liquid_solver.cpp)
target_include_directories(gifpy_liquid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(MSVC)
  target_compile_options(gifpy_liquid PRIVATE /W4 /permissive-)
else()
  target_compile_options(gifpy_liquid PRIVATE -Wall -Wextra -Wpedantic -Wconversion)
endif()

add_executable(gifpy_demo src/main.cpp)
target_link_libraries(gifpy_demo PRIVATE gifpy_liquid)

add_executable(gifpy_benchmark benchmarks/liquid_benchmark.cpp)
target_link_libraries(gifpy_benchmark PRIVATE gifpy_liquid)

enable_testing()
add_executable(gifpy_tests tests/liquid_solver_tests.cpp)
target_link_libraries(gifpy_tests PRIVATE gifpy_liquid)
add_test(NAME gifpy_liquid_tests COMMAND gifpy_tests)

add_test(NAME gifpy_liquid_benchmark_smoke COMMAND gifpy_benchmark --counts 64,256 --steps 3)
