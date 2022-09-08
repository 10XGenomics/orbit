cmake_minimum_required(VERSION 3.23)
project(STAR)

set(CMAKE_CXX_STANDARD 17)
add_subdirectory(source)
add_executable(STAR main.cpp)
target_link_libraries(STAR PUBLIC orbit)
target_include_directories(STAR PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/source"
        )