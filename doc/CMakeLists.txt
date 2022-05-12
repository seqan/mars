# ------------------------------------------------------------------------------------------------------------
# This is MaRs, Motif-based aligned RNA searcher.
# Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at https://github.com/seqan/mars.
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.8)

# Find doxygen.
find_package (Doxygen QUIET)

if (${DOXYGEN_FOUND})
    message (STATUS "Found Doxygen: ${DOXYGEN_EXECUTABLE}")

    # Configure doxygen options.
    set (MARS_DOXYFILE_IN ${CMAKE_SOURCE_DIR}/doc/doxygen_cfg)
    set (MARS_DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/doc")
    set (MARS_DOXYGEN_EXCLUDE_SYMBOLS "")
    set (MARS_DOXYGEN_PREDEFINED_NDEBUG "")
    set (MARS_DOXYGEN_ENABLED_SECTIONS "DEV")
    set (MARS_DOXYGEN_EXTRACT_PRIVATE "YES")

    configure_file (${MARS_DOXYFILE_IN} ${MARS_DOXYGEN_OUTPUT_DIR}/Doxyfile)
    add_custom_target (doc
                       COMMAND ${DOXYGEN_EXECUTABLE}
                       WORKING_DIRECTORY ${MARS_DOXYGEN_OUTPUT_DIR}
                       COMMENT "Generating (developer) API documentation with Doxygen."
                       VERBATIM)
    message (STATUS "${FontBold}You can run `make doc` to build api documentation.${FontReset}")
else ()
    message (STATUS "Doxygen not found.")
endif ()
