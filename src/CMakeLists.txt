cmake_minimum_required (VERSION 3.8)

# io
add_library (io STATIC multiple_alignment.cpp index_io.cpp)
target_link_libraries (io PUBLIC seqan3::seqan3)
target_include_directories (io PUBLIC .)

# structure
add_library (structure STATIC structure.cpp)
target_link_libraries (structure PUBLIC IPknot seqan3::seqan3)
target_include_directories (structure PUBLIC .)

# motif store
add_library (motif STATIC motif.cpp index.cpp search.cpp)
target_link_libraries (motif PUBLIC seqan3::seqan3 pthread)
target_include_directories (motif PUBLIC .)

# The mars executable consists of main.cpp and the linked object library.
add_executable ("${PROJECT_NAME}" main.cpp settings.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3 motif io structure pthread)
target_include_directories ("${PROJECT_NAME}" PUBLIC .)
