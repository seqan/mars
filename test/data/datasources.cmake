cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/
declare_datasource (FILE tRNA.aln
                    URL_HASH MD5=3508f9246b00552868f06bf26b3ef688
                    URL ${CMAKE_SOURCE_DIR}/test/data/tRNA.aln)

declare_datasource (FILE genome.fa
                    URL_HASH MD5=0b008b60b471145533f90ca1959ecf16
                    URL ${CMAKE_SOURCE_DIR}/test/data/genome.fa)

declare_datasource (FILE genome2.fa.marsindex
                    URL_HASH MD5=71457b32526a25fb0228a9fe152ee00b
                    URL ${CMAKE_SOURCE_DIR}/test/data/genome2.fa.marsindex)

declare_datasource (FILE genome3.fa.marsindex.gz
                    URL_HASH MD5=63cd20d19c9d30051203c6c049c4146d
                    URL ${CMAKE_SOURCE_DIR}/test/data/genome3.fa.marsindex.gz)

declare_datasource (FILE RF00005.fa.marsindex
                    URL_HASH MD5=26253dc665c78c6b80dc7a9024358bfd
                    URL ${CMAKE_SOURCE_DIR}/test/data/RF00005.fa.marsindex)

