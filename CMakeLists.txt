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
include (FindPkgConfig)
pkg_check_modules (ViennaRNA RNAlib2)

find_library (GLPK glpk)
find_library (MPFR mpfr)

find_package (OpenMP)
find_package (Boost)
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

if (NOT ViennaRNA_FOUND)
    include (ExternalProject)
    ExternalProject_Add (
            "ViennaRNA"
            PREFIX "${CMAKE_CURRENT_BINARY_DIR}/ext/ViennaRNA"
            LOG_DOWNLOAD 1
            DOWNLOAD_NO_PROGRESS 1
            TIMEOUT 240
            URL https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.4.17/ViennaRNA-2.4.17.tar.gz
            URL_HASH MD5=530b5ce7caf0e7cf58535cad22265215
            LOG_CONFIGURE 1
            USES_TERMINAL_CONFIGURE 1
            CONFIGURE_COMMAND <SOURCE_DIR>/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
                --prefix=<INSTALL_DIR> --without-swig --without-kinfold
                --without-forester --without-rnalocmin --without-doc --without-tutorial
            LOG_BUILD 1
            USES_TERMINAL_BUILD 1
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make -j
            LOG_INSTALL 1
    )
    ExternalProject_Get_Property (ViennaRNA INSTALL_DIR)
endif()

set (DEPENDENCIES_FOUND TRUE)

# Provided dependencies.
add_library (Contrafold STATIC lib/ipknot/contrafold/SStruct.cpp lib/ipknot/contrafold/Utilities.cpp)
target_include_directories (Contrafold SYSTEM PUBLIC lib/ipknot/contrafold)

add_library (Nupack STATIC lib/ipknot/nupack/nupack.cpp)
target_include_directories (Nupack SYSTEM PUBLIC lib/ipknot/nupack)
if (Boost_FOUND)
    target_include_directories (Nupack SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
else()
    message (WARNING "Boost not found, but it is required for building.")
    set (DEPENDENCIES_FOUND FALSE)
endif()

add_library (IPknot STATIC lib/ipknot/aln.cpp lib/ipknot/fold.cpp lib/ipknot/ip.cpp lib/ipknot/ipknot.cpp)
target_include_directories (IPknot SYSTEM PUBLIC lib/ipknot)
target_link_libraries (IPknot PUBLIC Contrafold Nupack)
target_compile_definitions (IPknot PRIVATE "-DHAVE_LIBRNA" PRIVATE "-DHAVE_VIENNA20")
target_compile_options(IPknot PRIVATE "-w")
if (ViennaRNA_FOUND)
    target_include_directories (IPknot SYSTEM PUBLIC ${ViennaRNA_INCLUDE_DIRS})
    target_link_libraries (IPknot PUBLIC ${ViennaRNA_LDFLAGS})
else ()
    add_dependencies (IPknot ViennaRNA)
    target_include_directories (IPknot SYSTEM PUBLIC ${INSTALL_DIR}/include)
    target_link_libraries (IPknot PUBLIC ${INSTALL_DIR}/lib/libRNA.a)
endif ()

# GLPK is needed for IPknot.
if (GLPK)
    target_compile_definitions (IPknot PRIVATE "-DWITH_GLPK")
    target_link_libraries (IPknot PUBLIC ${GLPK})
else()
    set (DEPENDENCIES_FOUND FALSE)
endif ()

# MPFR is needed for IPknot.
if (MPFR)
    target_link_libraries (IPknot PUBLIC ${MPFR})
else()
    set (DEPENDENCIES_FOUND FALSE)
endif ()

# OpenMP is needed for IPknot.
if (OpenMP_CXX_FOUND)
    target_link_libraries (IPknot PUBLIC OpenMP::OpenMP_CXX)
else()
    message (STATUS "OpenMP_CXX not found.")
    set (DEPENDENCIES_FOUND FALSE)
endif ()

add_subdirectory (src)

if (DEPENDENCIES_FOUND)
    message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")
else()
    message (STATUS "At least one dependency is missing.")
endif ()

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
message (STATUS "The current build type is ${CMAKE_BUILD_TYPE}.")
