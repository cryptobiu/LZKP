cmake_minimum_required(VERSION 3.5)
project(LZKP)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")

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

add_executable(LZKP
        test/unit_tests.cpp
        prover.cpp
        prover.h
        settings.h
        seedtree.cpp
        verifier.cpp
        verifier.h
        seedtree.h settings.h)

target_link_libraries(LZKP libcryptoTools.a libmiracl.a gsl ntl)