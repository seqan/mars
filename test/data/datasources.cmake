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
                    URL_HASH MD5=385e699d95b7e6ee6b5f33249b87dab9
                    URL ${CMAKE_SOURCE_DIR}/test/data/genome2.fa.marsindex)

declare_datasource (FILE genome3.fa.marsindex.gz
                    URL_HASH MD5=1e990d55f067fc95debd748d32953242
                    URL ${CMAKE_SOURCE_DIR}/test/data/genome3.fa.marsindex.gz)

declare_datasource (FILE RF0005.fa.marsindex
                    URL_HASH MD5=e0a6cb45c2ec8b50fd97bbe8aa56ddb5
                    URL ${CMAKE_SOURCE_DIR}/test/data/RF0005.fa.marsindex)

