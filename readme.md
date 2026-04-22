# F-Harmonic Forms on Symplectic Lie Groups - Computational Library

A C++ helper library for computing F-harmonic forms on 6-dimensional symplectic Lie algebras.

## Overview

This project implements algorithms to characterize **F-harmonic 3-forms** and analyze the **Type IIA flow** on symplectic Lie groups. It computes Lie algebra cohomology, symplectic cohomology, and polynomial invariants ($K$, $F$, $Q$ operators) for all 23 known 6-dimensional symplectic Lie algebras.

### Key Concepts

- **F-harmonic forms**: 3-forms $\phi$ satisfying both $d\phi = 0$ and $dF(\phi) = 0$
- **Symplectic Lie algebras**: Left-invariant symplectic structures on Lie groups
- **Primitive forms**: Forms orthogonal to the symplectic form under wedge product
- **Cohomology groups**: De Rham, primitive, and symplectic cohomologies
- **Type IIA flow**: Evolution equation $\frac{\partial \phi}{\partial t} = d\Lambda_\omega dF(\phi)$

## Mathematical Content

For detailed mathematical background and theory, see **[docs/MATHEMATICAL_CONTENT.md](docs/MATHEMATICAL_CONTENT.md)** which includes:

- Lie algebra foundations and structure constants
- Exterior derivative and Lie algebra cohomology
- Symplectic structures and primitive forms  
- Lefschetz decomposition and isomorphism
- Explicit K, F, Q operator definitions
- The main research questions and solution approaches

## Project Structure


### Directory Layout

```
.
├── include/Lie-Alg/           # Header files
│   ├── DifferentialForm.h     # Core exterior algebra operations
│   ├── Polynomials.h          # K, F, Q polynomial operators  
│   └── PairUtils.h            # Input parsing utilities
├── src/                       # Implementation
│   ├── main.cpp              # Process all 23 symplectic algebras
│   ├── DifferentialForm.cpp  # Wedge product, exterior derivative
│   └── Polynomials.cpp       # Explicit polynomial computations
├── test/                      # Unit tests
│   ├── test_differential_form.cpp  # DifferentialForm tests (15 tests)
│   ├── test_lie_algebra.cpp        # LieAlgebra tests (9 tests)
│   ├── test_framework.h            # Lightweight test framework
│   └── README.md                   # Testing documentation
├── docs/                      # Documentation
│   ├── MATHEMATICAL_CONTENT.md     # Detailed mathematical theory
│   └── symplectic-lie-groups.pdf   # Reference paper
├── res/                       # Input data
│   └── input.txt             # 23 symplectic Lie algebras
├── CMakeLists.txt            # Build configuration
├── Makefile                  # Convenience build target
└── output.md                 # Generated results
```

## Building the Project

### Prerequisites

- C++20 compiler
- GiNaC (symbolic computation library)
- CMake 3.29+

### Build Steps

```bash
cd /home/alfredo/Projects/Lie-Algebra
mkdir -p build
cd build
cmake ..
make
```

### Running the Program

```bash
./build/LIE-ALG
```

This processes all 23 symplectic Lie algebras from `res/input.txt` and generates results in `output.md`.

## Testing

This project includes a comprehensive unit test suite with **24 tests** covering:

- **DifferentialForm class**: constructors, wedge products, exterior derivatives, sign handling
- **LieAlgebra class**: structure constants, cohomology properties, nilpotence ($d^2 = 0$)
- **Mathematical properties**: Leibniz rule, additivity, Lefschetz decomposition

### Running Tests

```bash
cd build

# Method 1: Direct execution
./test/RunTests

# Method 2: Using CTest
ctest

# For detailed test documentation
cat ../test/README.md
```

## Implementation Details

### Core Components

**DifferentialForm**
- Represents exterior algebra elements as sparse polynomials
- Implements wedge product with automatic index deduplication
- Computes exterior derivatives using structure constants
- Handles GiNaC symbolic coefficients for polynomial computations

**LieAlgebra**
- Stores structure constants for a Lie algebra
- Provides differential operators for all basis elements
- Supports both general and primitive form computations

**Polynomials**
- Implements the K, F, Q operators as GiNaC expressions
- Precomputes all 20 possible 3-form basis elements
- Maps primitive and non-primitive combinations

### Input Format

The `res/input.txt` file contains 23 symplectic Lie algebra definitions in the format:

```
<name> <structure constants line 1> <structure constants line 2>
```

Each line uses `<index><coefficient>` pairs defining $d(e^i)$ relationships.

### Output

The program generates detailed results in `output.md` for each Lie algebra:
- Structure constants verification
- Symplectic form validity checking  
- Cohomology computations
- F-harmonic form characterization

## Key Features

✅ **Symbolic Computation** - Uses GiNaC for exact polynomial algebra

✅ **Comprehensive Testing** - 24 unit tests with custom lightweight framework

✅ **All 23 Algebras** - Exhaustive analysis of known symplectic Lie algebras

✅ **Efficient Representation** - Sparse polynomial storage for exterior forms

✅ **Documented Theory** - Complete mathematical background and proofs

## Performance Considerations

- Symbolic computations with large polynomials may be intensive
- Caching is implemented for structure constant derivatives
- Exterior product uses hash-based deduplication

## Contributing

To extend the project:

1. Add new test cases in `test/test_*.cpp`
2. Implement additional operators in `src/Polynomials.cpp`
3. Update documentation in `docs/MATHEMATICAL_CONTENT.md`
4. Submit changes on a feature branch

## References

The mathematical framework is based on:
- Symplectic Lie group classifications from `docs/symplectic-lie-groups.pdf`
- Lie algebra cohomology theory
- Hodge theory for symplectic structures

## Author Notes

This implementation focuses on:
- **Q1: Existence and uniqueness of F-harmonic representatives** in fixed cohomology classes
- **Q2: Long-time behavior** of the Type IIA evolution flow

All 23 symplectic Lie algebras can be analyzed systematically to answer these questions