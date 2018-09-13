cmake_minimum_required(VERSION 3.5)
project(LZKP)

set(CMAKE_CXX_STANDARD 14)

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("-- Compiliing in debug mode")
    add_definitions(-DDEBUG)
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic -g -O3 -march=native -mbmi -mbmi2")

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
        test/unit_tests_field15bit.cpp
        test/unit_tests_ZpMersenneIntElement.cpp
        test/unit_tests_ZpMersenneLongElement.cpp
        block.h
        utils.h
        parameters.h
        seedtree.cpp seedtree.h
        fields/field_15_bit.cpp fields/field_15_bit.h
        fields/mersenne.cpp fields/mersenne.h
        cac_prover_logic.cpp cac_prover_logic.h
        cac_prover.cpp cac_prover.h
        cac_verifier_logic.cpp cac_verifier_logic.h
        cac_verifier.cpp cac_verifier.h
        sac_prover_logic.cpp sac_prover_logic.h
        sac_prover.cpp sac_prover.h
        sac_verifier_logic.cpp sac_verifier_logic.h
        sac_verifier.cpp sac_verifier.h
        )

target_link_libraries(LZKP_TESTS libcryptoTools.a libmiracl.a gsl pthread gmp)

add_executable(LZKP
        main.cpp
        utils.h
        block.h
        parameters.h
        seedtree.cpp seedtree.h
        fields/field_15_bit.cpp fields/field_15_bit.h
        fields/mersenne.cpp fields/mersenne.h
        cac_prover_logic.cpp cac_prover_logic.h
        cac_prover.cpp cac_prover.h
        cac_verifier_logic.cpp cac_verifier_logic.h
        cac_verifier.cpp cac_verifier.h
        sac_prover_logic.cpp sac_prover_logic.h
        sac_prover.cpp sac_prover.h
        sac_verifier_logic.cpp sac_verifier_logic.h
        sac_verifier.cpp sac_verifier.h
        factory.cpp factory.h
        party.cpp party.h
        prover_party.cpp prover_party.h
        verifier_party.cpp verifier_party.h
        cac_prover_party.cpp cac_prover_party.h
        cac_verifier_party.cpp cac_verifier_party.h
        sac_prover_party.cpp sac_prover_party.h
        sac_verifier_party.cpp sac_verifier_party.h
        )

target_link_libraries(LZKP libcryptoTools.a libmiracl.a gsl pthread gmp boost_program_options)
