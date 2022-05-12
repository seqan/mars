# ------------------------------------------------------------------------------------------------------------
# This is MaRs, Motif-based aligned RNA searcher.
# Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at https://github.com/seqan/mars.
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.8)

add_api_test (bi_alphabet_test.cpp)

add_api_test (index_test.cpp)
target_use_datasources (index_test FILES
        genome.fa
        genome2.fa.marsindex
        genome3.fa.marsindex.gz
        RF00005.fa.marsindex
)

add_api_test (input_test.cpp)
target_use_datasources (input_test FILES tRNA.aln)
target_use_datasources (input_test FILES SSU_rRNA_5.sth)

add_api_test (motif_test.cpp)

add_api_test (profile_test.cpp)
