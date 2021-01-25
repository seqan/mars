cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/
declare_datasource (FILE tRNA.aln
                    URL ${CMAKE_SOURCE_DIR}/test/data/tRNA.aln)

declare_datasource (FILE genome.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/genome.fa)

declare_datasource (FILE genome2.fa.marsindex
                    URL ${CMAKE_SOURCE_DIR}/test/data/genome2.fa.marsindex)

declare_datasource (FILE RF0005.fa.marsindex
                    URL ${CMAKE_SOURCE_DIR}/test/data/RF0005.fa.marsindex)

