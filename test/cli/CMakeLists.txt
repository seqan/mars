cmake_minimum_required (VERSION 3.8)

add_cli_test (argument_parser_test.cpp)
target_use_datasources (argument_parser_test FILES tRNA.aln)
