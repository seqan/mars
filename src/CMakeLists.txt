cmake_minimum_required (VERSION 3.8)

# io
add_library (io STATIC input_output.cpp)
target_link_libraries (io PUBLIC seqan3::seqan3)
target_include_directories (io PUBLIC ../include)

# structure
add_library (structure STATIC structure.cpp)
target_link_libraries (structure PUBLIC IPknot seqan3::seqan3)
target_include_directories (structure PUBLIC ../include)

# motif store
add_library (motif STATIC motif.cpp)
target_link_libraries (motif PUBLIC seqan3::seqan3)
target_include_directories (motif PUBLIC ../include)

# index and search
add_library (index STATIC index.cpp)
target_link_libraries (index PUBLIC seqan3::seqan3)
target_include_directories (index PUBLIC ../include)

# The mars executable consists of main.cpp and the linked object library.
add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC motif io structure index)
