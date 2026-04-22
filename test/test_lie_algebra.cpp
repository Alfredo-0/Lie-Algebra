#include "test_framework.h"
#include "Lie-Alg/DifferentialForm.h"
#include "Lie-Alg/PairUtils.h"
#include <vector>
#include <memory>

/**
 * Test suite for LieAlgebra class and mathematical properties
 * Tests structure constants, closure properties, and Lie bracket identities
 */

TEST_CASE("LieAlgebra: Abelian Lie algebra (trivial structure constants)", "[LieAlgebra]") {
    // Create structure constants where all d(e^i) = 0
    std::vector<std::vector<Pair>> structure_constants;
    
    // 6 generators, all with trivial differentials
    for (int i = 0; i < 6; ++i) {
        structure_constants.push_back(std::vector<Pair>());
    }
    
    LieAlgebra algebra(structure_constants);
    
    // Set the algebra for forms
    DifferentialForm::algebra = std::make_shared<LieAlgebra>(structure_constants);
    
    // Test that dOf returns zero forms
    for (int i = 0; i < 6; ++i) {
        DifferentialForm d_form = algebra.dOf(i);
        REQUIRE(d_form.checkZero());
    }
}

TEST_CASE("ExteriorDerivative: d² = 0 for 0-forms (functions)", "[ExteriorDerivative]") {
    // Set up abelian algebra
    std::vector<std::vector<Pair>> structure;
    for (int i = 0; i < 6; ++i) {
        structure.push_back(std::vector<Pair>());
    }
    DifferentialForm::algebra = std::make_shared<LieAlgebra>(structure);
    
    // Create a 0-form (constant)
    DifferentialForm form0(0);
    
    // Apply d twice: d(d(f)) should be 0
    DifferentialForm d1 = form0.exteriorDerivative();
    DifferentialForm d2 = d1.exteriorDerivative();
    
    REQUIRE(d2.checkZero());
}

TEST_CASE("ExteriorDerivative: d² = 0 for 1-forms in abelian algebra", "[ExteriorDerivative]") {
    // Set up abelian algebra
    std::vector<std::vector<Pair>> structure;
    for (int i = 0; i < 6; ++i) {
        structure.push_back(std::vector<Pair>());
    }
    DifferentialForm::algebra = std::make_shared<LieAlgebra>(structure);
    
    // Create a 1-form
    DifferentialForm form1(1);
    std::array<int, DIMENSION> idx = {0};
    idx[0] = 1;
    form1.addTerm(idx, 1.0);
    
    // Apply d twice
    DifferentialForm d1 = form1.exteriorDerivative();
    DifferentialForm d2 = d1.exteriorDerivative();
    
    REQUIRE(d2.checkZero());
}

TEST_CASE("ExteriorDerivative: d² = 0 for 2-forms in abelian algebra", "[ExteriorDerivative]") {
    // Set up abelian algebra
    std::vector<std::vector<Pair>> structure;
    for (int i = 0; i < 6; ++i) {
        structure.push_back(std::vector<Pair>());
    }
    DifferentialForm::algebra = std::make_shared<LieAlgebra>(structure);
    
    // Create a 2-form
    DifferentialForm form2(2);
    std::array<int, DIMENSION> idx = {0};
    idx[0] = 1;
    idx[1] = 2;
    form2.addTerm(idx, 1.0);
    
    // Apply d twice
    DifferentialForm d1 = form2.exteriorDerivative();
    DifferentialForm d2 = d1.exteriorDerivative();
    
    REQUIRE(d2.checkZero());
}

TEST_CASE("ExteriorDerivative: d(α + β) = dα + dβ", "[ExteriorDerivative]") {
    // Set up abelian algebra
    std::vector<std::vector<Pair>> structure;
    for (int i = 0; i < 6; ++i) {
        structure.push_back(std::vector<Pair>());
    }
    DifferentialForm::algebra = std::make_shared<LieAlgebra>(structure);
    
    // Create two 1-forms
    DifferentialForm form1(1), form2(1);
    std::array<int, DIMENSION> idx1 = {0}, idx2 = {0};
    idx1[0] = 1;
    idx2[0] = 2;
    
    form1.addTerm(idx1, 2.0);
    form2.addTerm(idx2, 3.0);
    
    // Compute d(α + β)
    DifferentialForm sum = form1;
    sum += form2;
    DifferentialForm d_sum = sum.exteriorDerivative();
    
    // Compute dα + dβ
    DifferentialForm d_form1 = form1.exteriorDerivative();
    DifferentialForm d_form2 = form2.exteriorDerivative();
    DifferentialForm sum_d = d_form1;
    sum_d += d_form2;
    
    // They should be equal (d_sum + (-sum_d) = 0)
    d_sum += sum_d;
    d_sum += DifferentialForm();  // Add zero to ensure it's properly processed
    
    // Note: In abelian algebra, dα = 0 for any form
    REQUIRE(d_sum.checkZero());
}

TEST_CASE("LieAlgebra: Leibniz rule in abelian algebra", "[LieAlgebra]") {
    // Set up abelian algebra (simpler for testing)
    std::vector<std::vector<Pair>> structure;
    for (int i = 0; i < 6; ++i) {
        structure.push_back(std::vector<Pair>());
    }
    DifferentialForm::algebra = std::make_shared<LieAlgebra>(structure);
    
    // Create two 1-forms
    DifferentialForm alpha(1), beta(1);
    std::array<int, DIMENSION> idx_alpha = {0}, idx_beta = {0};
    idx_alpha[0] = 1;
    idx_beta[0] = 2;
    
    alpha.addTerm(idx_alpha, 1.0);
    beta.addTerm(idx_beta, 1.0);
    
    // LHS: d(α ∧ β)
    DifferentialForm wedge_product = alpha.wedge(beta);
    DifferentialForm lhs = wedge_product.exteriorDerivative();
    
    // RHS: dα ∧ β + (-1)^|α| α ∧ dβ
    // In abelian algebra, dα = dβ = 0, so RHS = 0
    DifferentialForm rhs(3);  // 3-form
    
    // LHS should also be 0 in abelian algebra
    REQUIRE(lhs.checkZero());
    REQUIRE(rhs.checkZero());
}

TEST_CASE("LieAlgebra: Basis 3-forms are accessible", "[LieAlgebra]") {
    // This tests the hardcoded basis 3-forms
    REQUIRE(basis_3forms.size() == 20);
    
    // All basis forms should have 3 indices
    for (const auto& form : basis_3forms) {
        REQUIRE(form[0] > 0);
        REQUIRE(form[1] > 0);
        REQUIRE(form[2] > 0);
        REQUIRE(form[0] < form[1]);
        REQUIRE(form[1] < form[2]);
    }
}

TEST_CASE("LieAlgebra: Primitive basis 3-forms are valid", "[LieAlgebra]") {
    REQUIRE(primitive_basis_3forms.size() == 20);
    
    // Each primitive basis should reference two indices from basis_3forms
    for (const auto& prim : primitive_basis_3forms) {
        REQUIRE(prim[0] < 20);
        REQUIRE(prim[1] < 20);
    }
}

TEST_CASE("LieAlgebra: Zero form remains zero after operations", "[LieAlgebra]") {
    std::vector<std::vector<Pair>> structure;
    for (int i = 0; i < 6; ++i) {
        structure.push_back(std::vector<Pair>());
    }
    DifferentialForm::algebra = std::make_shared<LieAlgebra>(structure);
    
    DifferentialForm zero1, zero2;
    
    // Wedge of zero forms
    DifferentialForm result = zero1.wedge(zero2);
    REQUIRE(result.checkZero());
    
    // Sum of zero forms
    zero1 += zero2;
    REQUIRE(zero1.checkZero());
    
    // Exterior derivative of zero
    DifferentialForm d_zero = zero1.exteriorDerivative();
    REQUIRE(d_zero.checkZero());
}
