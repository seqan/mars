cmake_minimum_required (VERSION 3.8)

# External dependencies
find_library(ViennaRNA libRNA.a)
find_library(GLPK glpk)

# Provided dependencies.
add_library(Contrafold STATIC ../lib/ipknot/contrafold/SStruct.cpp ../lib/ipknot/contrafold/Utilities.cpp)
target_include_directories (Contrafold PUBLIC ../lib/ipknot/contrafold)

add_library(Nupack STATIC ../lib/ipknot/nupack/nupack.cpp)
target_include_directories (Nupack PUBLIC ../lib/ipknot/nupack)

add_library(IPknot STATIC ../lib/ipknot/aln.cpp ../lib/ipknot/fold.cpp ../lib/ipknot/ip.cpp ../lib/ipknot/ipknot.cpp)
target_include_directories (IPknot PUBLIC ../lib/ipknot)
target_link_libraries (IPknot PUBLIC Contrafold Nupack)

if (ViennaRNA)
    message(STATUS "Found ViennaRNA: ${ViennaRNA}")
    target_compile_definitions(IPknot PRIVATE "-DHAVE_LIBRNA" PRIVATE "-DHAVE_VIENNA18" PRIVATE "-DHAVE_VIENNA20")
    target_link_libraries (IPknot PUBLIC ${ViennaRNA})
endif ()

if (GLPK)
    message(STATUS "Found GLPK: ${GLPK}")
    target_compile_definitions(IPknot PRIVATE "-DWITH_GLPK")
    target_link_libraries (IPknot PUBLIC ${GLPK})
endif ()

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" OBJECT motif_store.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3 IPknot)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC "${PROJECT_NAME}_lib")
