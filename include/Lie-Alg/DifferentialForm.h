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

struct Triple{
    int i;
    int j;
    int k;
};

struct Pairs{
    int i;
    int j;
};

const std::array<Triple, 20> basis_3forms = {{
    {1,2,3}, {1,2,4}, {1,2,5}, {1,2,6}, 
    {1,3,4}, {1,3,5}, {1,3,6}, {1,4,5},
    {1,4,6}, {1,5,6}, {2,3,4}, {2,3,5},
    {2,3,6}, {2,4,5}, {2,4,6}, {2,5,6}, 
    {3,4,5}, {3,4,6}, {3,5,6}, {4,5,6}
}};

const std::array<Pairs, 15> basis_2forms = {{
    {1,2}, {1,3}, {1,4}, {1,5}, {1,6},
    {2,3}, {2,4}, {2,5}, {2,6}, {3,4},
    {3,5}, {3,6}, {4,5}, {4,6}, {5,6}
}};

class LieAlgebra;

class DifferentialForm {
private:
    int degree;
public:
    std::map<std::array<int, DIMENSION>, double> terms;//change!

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

    DifferentialForm inverse() const;

};

class LieAlgebra {
private:
    std::array<DifferentialForm, DIMENSION> structureConstants;

public:
    
    LieAlgebra(std::vector<std::vector<Pair>> str);

    DifferentialForm dOf(int i) const {
        return structureConstants[i-1];
    }
};