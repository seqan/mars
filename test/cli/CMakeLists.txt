# ------------------------------------------------------------------------------------------------------------
# This is MaRs, Motif-based aligned RNA searcher.
# Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at https://github.com/seqan/mars.
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.8)

add_cli_test (argument_parser_test.cpp)
target_use_datasources (argument_parser_test FILES tRNA.aln)
