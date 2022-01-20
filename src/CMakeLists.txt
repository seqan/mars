cmake_minimum_required (VERSION 3.8)

# settings (needs thread_pool)
add_library(settings STATIC settings.cpp)
target_include_directories(settings PUBLIC ../lib/thread_pool)
target_link_libraries(settings PUBLIC seqan3::seqan3 pthread)

# structural multiple alignment (needs IPknot)
add_library (msa STATIC multiple_alignment.cpp)
target_link_libraries (msa PUBLIC IPknot seqan3::seqan3)

# index io (needs settings)
add_library (index STATIC index.cpp)
target_link_libraries (index PUBLIC seqan3::seqan3 settings)

# motif (needs msa and settings)
add_library (motif STATIC motif.cpp)
target_link_libraries (motif PUBLIC seqan3::seqan3 msa settings)

# search
add_library (search STATIC search.cpp)
target_link_libraries (search PUBLIC seqan3::seqan3 motif index settings)

# The mars executable consists of main.cpp and the linked object library.
add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3 settings search index msa motif)
