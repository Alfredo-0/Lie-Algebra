cmake_minimum_required(VERSION 3.29)

# Define test executable
add_executable(RunTests
    test_main.cpp
    test_differential_form.cpp
    test_lie_algebra.cpp
)

# Link libraries
target_link_libraries(RunTests PRIVATE ${GINAC_LIBRARIES})

# Include directories
target_include_directories(RunTests PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(RunTests PRIVATE ${GINAC_INCLUDE_DIRS})

# Compiler flags
target_compile_options(RunTests PRIVATE ${GINAC_CFLAGS_OTHER})

# Link against DifferentialForm and Polynomials implementations
target_sources(RunTests PRIVATE 
    ${PROJECT_SOURCE_DIR}/src/DifferentialForm.cpp
    ${PROJECT_SOURCE_DIR}/src/Polynomials.cpp
)

# Register tests
enable_testing()
add_test(NAME RunTests COMMAND RunTests)
