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

add_api_test (motif_test.cpp)

add_api_test (profile_test.cpp)
