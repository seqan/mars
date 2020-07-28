cmake_minimum_required (VERSION 3.8)

add_api_test (input_test.cpp)
target_use_datasources (input_test FILES tRNA.aln)
