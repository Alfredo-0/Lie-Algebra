cmake_minimum_required (VERSION 3.29)

project(LIE-ALG VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

file(GLOB project_sources
 "src/*.cpp"
)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(${PROJECT_NAME} ${project_sources})

if(EXISTS "${PROJECT_SOURCE_DIR}/test/CMakeLists.txt")
    enable_testing()
    add_subdirectory(test)
endif()