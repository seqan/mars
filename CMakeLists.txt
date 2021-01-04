cmake_minimum_required (VERSION 3.8)

## CUSTOMISE

# Define the application name and version.
project (mars VERSION 1.0.0)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

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
find_library(GLPK glpk)
find_library(MPFR mpfr)
find_library(ViennaRNA libRNA.a)
find_package(OpenMP)
find_package(Boost)
find_package(SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

# Provided dependencies.
add_library(Contrafold STATIC lib/ipknot/contrafold/SStruct.cpp lib/ipknot/contrafold/Utilities.cpp)
target_include_directories (Contrafold SYSTEM PUBLIC lib/ipknot/contrafold)

add_library(Nupack STATIC lib/ipknot/nupack/nupack.cpp)
target_include_directories (Nupack SYSTEM PUBLIC lib/ipknot/nupack)

add_library(IPknot STATIC lib/ipknot/aln.cpp lib/ipknot/fold.cpp lib/ipknot/ip.cpp lib/ipknot/ipknot.cpp)
target_include_directories (IPknot SYSTEM PUBLIC lib/ipknot)
target_link_libraries (IPknot PUBLIC Contrafold Nupack)

add_subdirectory (src)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC include)

# Boost is needed for Nupack and Contrafold.
if (Boost_FOUND)
    target_include_directories(Nupack PUBLIC ${Boost_INCLUDE_DIRS})
    target_include_directories(Contrafold PUBLIC ${Boost_INCLUDE_DIRS})
endif()

# ViennaRNA is needed for IPknot.
if (ViennaRNA)
    message(STATUS "Found ViennaRNA: ${ViennaRNA}")
    target_compile_definitions(IPknot PRIVATE "-DHAVE_LIBRNA" PRIVATE "-DHAVE_VIENNA20")
    target_link_libraries (IPknot PUBLIC ${ViennaRNA})
else()
    message(STATUS "ViennaRNA not found.")
endif ()

# MPFR is needed for IPknot if Vienna was compiled with MPFR.
if (MPFR)
    message(STATUS "Found MPFR: ${MPFR}")
    target_link_libraries (IPknot PUBLIC ${MPFR})
else()
    message(STATUS "MPFR not found.")
endif ()

# GLPK is needed for IPknot.
if (GLPK)
    message(STATUS "Found GLPK: ${GLPK}")
    target_compile_definitions(IPknot PRIVATE "-DWITH_GLPK")
    target_link_libraries (IPknot PUBLIC ${GLPK})
else()
    message(STATUS "GLPK not found.")
endif ()

# OpenMP is needed for IPknot.
if (OpenMP_CXX_FOUND)
    target_link_libraries(IPknot PUBLIC OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP not found.")
endif ()

message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
message(STATUS "The current build type is ${CMAKE_BUILD_TYPE}.")
