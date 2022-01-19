cmake_minimum_required (VERSION 3.8)

# settings (needs thread_pool)
add_library(settings STATIC settings.cpp)
target_include_directories(settings PUBLIC ../lib/thread_pool)
target_link_libraries(settings PUBLIC seqan3::seqan3 pthread)

# index io
add_library (io STATIC index_io.cpp)
target_link_libraries (io PUBLIC seqan3::seqan3)

# structural multiple alignment (needs IPknot)
add_library (msa STATIC multiple_alignment.cpp)
target_link_libraries (msa PUBLIC IPknot seqan3::seqan3)

# motif store
add_library (motif STATIC motif.cpp index.cpp search.cpp)
target_link_libraries (motif PUBLIC seqan3::seqan3 settings io msa)

# The mars executable consists of main.cpp and the linked object library.
add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3 settings io msa motif)
