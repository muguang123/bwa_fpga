cmake_minimum_required(VERSION 2.9)
project(bwa_fpga)

set(CMAKE_CXX_STANDARD 11)

include_directories(include include/junction include/turf)
link_directories(libs)
find_library(zlib z)

link_libraries(pthread libjunction.a libturf.a )

aux_source_directory(./ SRC_FILES)

#message("${SRC_FILES}")

add_executable(bwa_fpga ${SRC_FILES})


target_link_libraries(bwa_fpga ${zlib} m rt)
