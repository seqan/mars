cmake_minimum_required (VERSION 3.8)

add_api_test (bi_alphabet_test.cpp)

add_api_test (input_test.cpp)
target_use_datasources (input_test FILES tRNA.aln)

add_api_test (motif_store_test.cpp)

add_api_test (profile_test.cpp)
