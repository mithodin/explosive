cmake_minimum_required(VERSION 2.8.1)
project(explosive)

#set(CMAKE_BUILD_TYPE Debug)

include_directories("/opt/intel/include")
link_directories("/opt/intel/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")
set(CMAKE_C_FLAGS "-O3 -march=core-avx2 -std=gnu99")
#set(CMAKE_C_FLAGS "-O0 -march=core2 -std=gnu99")
add_definitions(-DDSFMT_MEXP=19937 -DHAVE_SSE2)
add_executable(explosive dSFMT/dSFMT.c explosive.c geometry.c colloid.c monte_carlo.c hdf5_logging.c logger.c substrate.c clusters.c helpers.c)
target_link_libraries(explosive hdf5 hdf5_hl pthread)

set_property( TARGET explosive APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )

cmake_policy(VERSION 2.8.1)
