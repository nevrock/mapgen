cmake_minimum_required(VERSION 3.0)
project(mapgen)

include_directories(${CMAKE_SOURCE_DIR}/include/)

add_library(mapgen SHARED src/stb_image.cpp src/mapgen.cpp)

# Add the executable
add_executable(mapgen_exe src/stb_image.cpp src/main.cpp) 

# (Optional) Link the executable to the shared library
# This is only needed if your executable uses functions from the library
# target_link_libraries(mapgen_exe mapgen) 