cmake_minimum_required(VERSION 3.5)
project(LZKP)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic -g -O3 -march=native")

include_directories(.)

link_directories(/usr/local/lib /usr/lib)

find_package(Boost 1.62.0 REQUIRED
        #[version] [EXACT]      # Minimum or EXACT version e.g. 1.36.0
        #[REQUIRED]             # Fail with error if Boost is not found
        #[COMPONENTS <libs>...] # Boost libraries by their canonical name
        )

IF(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
    LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
ENDIF(Boost_FOUND)

add_executable(LZKP_TESTS
        test/unit_tests.cpp
        settings.h
        prover_wrapper.cpp
        prover_wrapper.h
        prover.cpp
        prover.h
        verifier_wrapper.cpp
        verifier_wrapper.h
        verifier.cpp
        verifier.h
        seedtree.cpp
        seedtree.h
        Mersenne.cpp
        Mersenne.h
        Ring31.cpp
        Ring31.h
        RingNTL.cpp
        RingNTL.h)

target_link_libraries(LZKP_TESTS libcryptoTools.a libmiracl.a gsl ntl pthread gmp)

add_executable(LZKP
        protocol.cpp
        settings.h
        prover_wrapper.cpp
        prover_wrapper.h
        prover.cpp
        prover.h
        verifier_wrapper.cpp
        verifier_wrapper.h
        verifier.cpp
        verifier.h
        seedtree.cpp
        seedtree.h
        Mersenne.cpp
        Mersenne.h)

target_link_libraries(LZKP libcryptoTools.a libmiracl.a gsl ntl pthread gmp)