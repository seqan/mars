# ------------------------------------------------------------------------------------------------------------
# This is MaRs, Motif-based aligned RNA searcher.
# Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at https://github.com/seqan/mars.
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.13)

# Create an object library that is shared between main and tests.
add_library(lib${PROJECT_NAME} OBJECT
        index.cpp
        location.cpp
        motif.cpp
        multiple_alignment.cpp
        search.cpp
        settings.cpp
)
target_include_directories(lib${PROJECT_NAME} PUBLIC ../lib/thread_pool)
target_link_libraries (lib${PROJECT_NAME} PUBLIC seqan3::seqan3 IPknot pthread)

# Create the main executable.
add_executable (${PROJECT_NAME} main.cpp $<TARGET_OBJECTS:lib${PROJECT_NAME}>)
target_include_directories(${PROJECT_NAME} PUBLIC ../lib/thread_pool)
target_link_libraries (${PROJECT_NAME} PUBLIC seqan3::seqan3 IPknot pthread)
