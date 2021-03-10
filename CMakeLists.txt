cmake_minimum_required (VERSION 3.8)

## CUSTOMISE

# Define the application name and version.
project (mars VERSION 1.0.0)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# External dependencies
# Alternative for PkgConfig: https://dominikberner.ch/cmake-find-library/
include(FindPkgConfig)
pkg_check_modules(ViennaRNA RNAlib2)

find_library(GLPK glpk)
find_library(MPFR mpfr)

find_package(OpenMP)
find_package(Boost)
find_package(SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

# Provided dependencies.
add_library(Contrafold STATIC lib/ipknot/contrafold/SStruct.cpp lib/ipknot/contrafold/Utilities.cpp)
target_include_directories (Contrafold SYSTEM PUBLIC lib/ipknot/contrafold)

add_library(Nupack STATIC lib/ipknot/nupack/nupack.cpp)
target_include_directories (Nupack SYSTEM PUBLIC lib/ipknot/nupack)
if (Boost_FOUND)
    target_include_directories (Nupack SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
endif()

add_library(IPknot STATIC lib/ipknot/aln.cpp lib/ipknot/fold.cpp lib/ipknot/ip.cpp lib/ipknot/ipknot.cpp)
target_include_directories (IPknot SYSTEM PUBLIC lib/ipknot)
target_link_libraries (IPknot PUBLIC Contrafold Nupack)
if (ViennaRNA_FOUND)
    target_include_directories (IPknot SYSTEM PUBLIC ${ViennaRNA_INCLUDE_DIRS})
    target_link_libraries (IPknot PUBLIC ${ViennaRNA_LDFLAGS} ${MPFR})
    target_compile_definitions(IPknot PRIVATE "-DHAVE_LIBRNA" PRIVATE "-DHAVE_VIENNA20")
else ()
    message(WARNING "ViennaRNA not found: Builds will fail.")
endif ()

# GLPK is needed for IPknot.
if (GLPK)
    target_compile_definitions(IPknot PRIVATE "-DWITH_GLPK")
    target_link_libraries (IPknot PUBLIC ${GLPK})
endif ()

# OpenMP is needed for IPknot.
if (OpenMP_CXX_FOUND)
    target_link_libraries(IPknot PUBLIC OpenMP::OpenMP_CXX)
    add_definitions(-DMARS_WITH_OPENMP)
else()
    message(STATUS "OpenMP_CXX not found.")
endif ()

add_subdirectory (src)

message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
message(STATUS "The current build type is ${CMAKE_BUILD_TYPE}.")
