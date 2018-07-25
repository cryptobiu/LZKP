cmake_minimum_required(VERSION 3.5)
project(LZKP)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(LZKP
        test/unit_tests.cpp
        prover.cpp
        prover.h)
