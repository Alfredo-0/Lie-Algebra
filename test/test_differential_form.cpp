#include "test_framework.h"
#include "Lie-Alg/DifferentialForm.h"
#include "Lie-Alg/PairUtils.h"
#include <array>

/**
 * Test suite for DifferentialForm class
 * Tests basic operations: addTerm, checkZero, wedge product, exterior derivative
 */

TEST_CASE("DifferentialForm: Default constructor creates zero form", "[DifferentialForm]") {
    DifferentialForm form;
    REQUIRE(form.checkZero());
}

TEST_CASE("DifferentialForm: Constructor with degree", "[DifferentialForm]") {
    DifferentialForm form1(1);
    DifferentialForm form2(2);
    DifferentialForm form3(3);
    REQUIRE(form1.checkZero());
    REQUIRE(form2.checkZero());
    REQUIRE(form3.checkZero());
}

TEST_CASE("DifferentialForm: Adding a single term creates non-zero form", "[DifferentialForm]") {
    DifferentialForm form;
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;  // e^1
    form.addTerm(indices, 2.0);
    
    REQUIRE(!form.checkZero());
}

TEST_CASE("DifferentialForm: Adding identical terms sums correctly", "[DifferentialForm]") {
    DifferentialForm form;
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;  // e^1
    
    form.addTerm(indices, 3.0);
    form.addTerm(indices, 2.0);
    
    REQUIRE(!form.checkZero());
}

TEST_CASE("DifferentialForm: Adding opposite coefficients cancels", "[DifferentialForm]") {
    DifferentialForm form;
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;  // e^1
    
    form.addTerm(indices, 5.0);
    form.addTerm(indices, -5.0);
    
    REQUIRE(form.checkZero());
}

TEST_CASE("DifferentialForm: Sign handling for basis 1-forms", "[DifferentialForm]") {
    DifferentialForm form1, form2;
    std::array<int, DIMENSION> indices1 = {0}, indices2 = {0};
    
    indices1[0] = 1; indices1[1] = 2;  // e^{1,2}
    indices2[0] = 2; indices2[1] = 1;  // e^{2,1}
    
    form1.addTerm(indices1, 1.0);
    form2.addTerm(indices2, 1.0);
    
    form1 += form2;
    REQUIRE(form1.checkZero());
}

TEST_CASE("DifferentialForm: Wedge product with zero form", "[DifferentialForm]") {
    DifferentialForm form1(1), form2;
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;  // e^1
    form1.addTerm(indices, 1.0);
    
    DifferentialForm result = form1.wedge(form2);
    REQUIRE(result.checkZero());
}

TEST_CASE("DifferentialForm: Wedge product of 1-forms", "[DifferentialForm]") {
    DifferentialForm form1(1), form2(1);
    
    std::array<int, DIMENSION> indices1 = {0}, indices2 = {0};
    indices1[0] = 1;  // e^1
    indices2[0] = 2;  // e^2
    
    form1.addTerm(indices1, 1.0);
    form2.addTerm(indices2, 1.0);
    
    DifferentialForm result = form1.wedge(form2);
    REQUIRE(!result.checkZero());
}

TEST_CASE("DifferentialForm: Repeated index in wedge product gives zero", "[DifferentialForm]") {
    DifferentialForm form1(1), form2(1);
    
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;  // e^1
    
    form1.addTerm(indices, 1.0);
    form2.addTerm(indices, 1.0);
    
    DifferentialForm result = form1.wedge(form2);
    REQUIRE(result.checkZero());
}

TEST_CASE("DifferentialForm: Anticommutativity of wedge product", "[DifferentialForm]") {
    DifferentialForm form1(1), form2(1);
    
    std::array<int, DIMENSION> idx1 = {0}, idx2 = {0};
    idx1[0] = 1;  // e^1
    idx2[0] = 2;  // e^2
    
    form1.addTerm(idx1, 1.0);
    form2.addTerm(idx2, 1.0);
    
    DifferentialForm w1 = form1.wedge(form2);
    DifferentialForm w2 = form2.wedge(form1);
    
    // w1 + w2 should be zero (anticommutativity: ω ∧ η = -η ∧ ω)
    w1 += w2;
    REQUIRE(w1.checkZero());
}

TEST_CASE("DifferentialForm: Forms with GiNaC symbols", "[DifferentialForm]") {
    initialize_symbols();
    
    DifferentialForm form;
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;  // e^1
    
    // Add a term with a symbolic coefficient
    form.addTerm(indices, A);
    
    REQUIRE(!form.checkZero());
}

TEST_CASE("DifferentialForm: operator+= adding forms", "[DifferentialForm]") {
    DifferentialForm form1, form2;
    std::array<int, DIMENSION> idx = {0};
    idx[0] = 1;
    
    form1.addTerm(idx, 2.0);
    form2.addTerm(idx, 3.0);
    
    form1 += form2;
    REQUIRE(!form1.checkZero());
}

TEST_CASE("DifferentialForm: operator+= adding and canceling forms", "[DifferentialForm]") {
    DifferentialForm form1, form2;
    std::array<int, DIMENSION> idx = {0};
    idx[0] = 1;
    
    form1.addTerm(idx, 5.0);
    form2.addTerm(idx, -5.0);
    
    form1 += form2;
    REQUIRE(form1.checkZero());
}

TEST_CASE("DifferentialForm: 1-form degree is 1", "[DifferentialForm]") {
    DifferentialForm form;
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;
    form.addTerm(indices, 1.0);
    
    // Degree should be set to 1
    // Note: degree is private, so we test through behavior
    DifferentialForm wedged = form.wedge(form);
    REQUIRE(wedged.checkZero());  // 1-form ∧ 1-form with same index = 0
}

TEST_CASE("DifferentialForm: 2-form degree is 2", "[DifferentialForm]") {
    DifferentialForm form;
    std::array<int, DIMENSION> indices = {0};
    indices[0] = 1;
    indices[1] = 2;
    form.addTerm(indices, 1.0);
    
    DifferentialForm form1;
    std::array<int, DIMENSION> idx1 = {0};
    idx1[0] = 3;
    form1.addTerm(idx1, 1.0);
    
    DifferentialForm result = form.wedge(form1);
    REQUIRE(!result.checkZero());  // 2-form ∧ 1-form = 3-form
}
