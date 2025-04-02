cmake_minimum_required (VERSION 3.29)

project(LIE-ALG VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

file(GLOB project_sources
 "src/*.cpp"
)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GINAC REQUIRED ginac)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(${PROJECT_NAME} ${project_sources})

target_compile_definitions("${CMAKE_PROJECT_NAME}" PUBLIC RESOURCES_PATH="${CMAKE_CURRENT_SOURCE_DIR}/res/input.txt")
target_compile_definitions("${CMAKE_PROJECT_NAME}" PUBLIC DOCS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/output.md")

target_include_directories(${PROJECT_NAME} PRIVATE ${GINAC_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${GINAC_LIBRARIES})
target_compile_options(${PROJECT_NAME} PRIVATE ${GINAC_CFLAGS_OTHER})

if(EXISTS "${PROJECT_SOURCE_DIR}/test/CMakeLists.txt")
    enable_testing()
    add_subdirectory(test)
endif()