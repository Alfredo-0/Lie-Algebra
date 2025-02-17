//DifferentialForm.h

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <array>
#include <map>

#include "Lie-Alg/PairUtils.h"

const int DIMENSION = 6;


class LieAlgebra;

class DifferentialForm {
private:
    std::map<std::array<int, DIMENSION>, double> terms;
    int degree;
public:

    inline static std::shared_ptr<LieAlgebra> algebra = nullptr;

    DifferentialForm() {}
    DifferentialForm(int d) : degree(d) { }

    ~DifferentialForm() {}

    void addTerm(const std::array<int, DIMENSION>& indices, double coeff);

    void print() const;

    std::string toLaTeX() const;

    bool checkZero() const;

    DifferentialForm wedge(const DifferentialForm& other) const;
    
    DifferentialForm exteriorDerivative() const;

    DifferentialForm& operator +=(const DifferentialForm& other);

    DifferentialForm interiorProduct(const DifferentialForm& other) const;

};

class LieAlgebra {
private:
    std::array<DifferentialForm, DIMENSION> structureConstants;

public:
    
    LieAlgebra(std::vector<Pair> str);

    DifferentialForm dOf(int i) const {
        return structureConstants[i-1];
    }
};