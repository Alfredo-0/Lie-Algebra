# Test Suite

This directory contains the unit tests for the Lie Algebra computation project.

## Building Tests

The tests use **Catch2** as the testing framework (header-only, automatically downloaded via CMake).

```bash
cd /home/alfredo/Projects/Lie-Algebra
mkdir -p build
cd build
cmake ..
make
```

## Running Tests

After building:

```bash
# Run all tests
ctest

# Or run the test executable directly for more verbose output
./RunTests

# Run tests with a specific tag
./RunTests "[DifferentialForm]"

# Run with verbose output
./RunTests -v
```

## Test Coverage

### test_differential_form.cpp
Tests for the `DifferentialForm` class:
- **Initialization**: Constructor behavior and zero forms
- **addTerm()**: Coefficient handling, sign management, term cancellation
- **checkZero()**: Zero detection
- **wedge()**: Wedge product rules, anticommutativity
- **Symbolic Integration**: GiNaC symbol handling
- **operator+=**: Form addition
- **Degree Tracking**: Automatic degree inference

### test_lie_algebra.cpp
Tests for `LieAlgebra` class and mathematical properties:
- **Lie Algebra Construction**: Structure constant initialization
- **Exterior Derivative Properties**: 
  - $d^2 = 0$ (nilpotence)
  - Additivity: $d(\alpha + \beta) = d\alpha + d\beta$
- **Leibniz Rule**: $d(\alpha \wedge \beta) = d\alpha \wedge \beta + (-1)^{|\alpha|} \alpha \wedge d\beta$
- **Primitive Forms**: Basis 3-form validation
- **Zero Properties**: Behavior under operations

## Key Mathematical Properties Tested

1. **Nilpotence**: $d^2 = 0$ for all degrees
2. **Linearity**: $d(c \cdot \alpha) = c \cdot d(\alpha)$
3. **Additivity**: $d(\alpha + \beta) = d\alpha + d\beta$
4. **Anticommutativity**: $\alpha \wedge \beta = -\beta \wedge \alpha$
5. **Repeated Index Rule**: $\alpha_i \wedge \alpha_i = 0$
6. **Sign Handling**: Proper tracking of signs in basis sorting

## Adding New Tests

To add new tests:

1. Add test cases to the appropriate file or create a new test file
2. Include `<catch2/catch_test_macros.hpp>`
3. Use the `TEST_CASE` macro with clear test names and sections
4. Rebuild with `make` from the build directory

Example:
```cpp
TEST_CASE("MyTest: Description", "[Tag]") {
    SECTION("Specific behavior") {
        // Test code here
        REQUIRE(condition);
    }
}
```

## Future Test Coverage

Planned tests:
- Interior product operations
- Inverse operations for symplectic forms
- LaTeX output formatting
- F-harmonic form computations
- Full symplectic Lie algebra examples
- Performance benchmarks
