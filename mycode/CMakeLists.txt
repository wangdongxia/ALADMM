cmake_minimum_required(VERSION 3.8)
project(RADADMM)

set(CMAKE_CXX_COMPILER "mpic++")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR}/include)
link_libraries(pthread)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_FILE)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/internal SRC_FILE)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/data SRC_FILE)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/math SRC_FILE)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/optimizer SRC_FILE)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/other SRC_FILE)

add_executable(admm ${SRC_FILE})
