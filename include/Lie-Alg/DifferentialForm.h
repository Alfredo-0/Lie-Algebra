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

const int DIMENSION = 12;

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

struct Comparator;
struct PairComparator;

class LieAlgebra;

class DifferentialForm {
    std::map<std::array<int, DIMENSION>, double> terms;
    int degree;
    bool degree_assigned;
    
    public:
    inline static std::shared_ptr<LieAlgebra> algebra = nullptr;

    DifferentialForm(int d) : degree(d), degree_assigned(true) { }
    
    DifferentialForm(const std::array<int, DIMENSION>& indices, double coeff)
    : degree(0), degree_assigned(false) {
        addTerm(indices, coeff);
    }
    
    DifferentialForm() : degree(0), degree_assigned(false) {}   

    ~DifferentialForm() {}

    void addTerm(const std::array<int, DIMENSION>& indices, double coeff);
    
    friend struct Comparator;
    
    std::string toLaTeX() const;
    
    bool checkZero() const;
    
    DifferentialForm wedge(const DifferentialForm& other) const;
    
    DifferentialForm exteriorDerivative() const;
    
    DifferentialForm& operator +=(const DifferentialForm& other);
    
    DifferentialForm interiorProduct(const DifferentialForm& other) const;
    
    DifferentialForm inverse() const;
    
    void print() const;
};

class LieAlgebra {
    std::array<DifferentialForm, DIMENSION> structureConstants;

public:
    LieAlgebra(std::vector<std::vector<Pair>> str);

    DifferentialForm dOf(int i) const {
        return structureConstants[i];
    }

};

struct Comparator {
    bool operator()(const DifferentialForm& a, const DifferentialForm b) const {
        auto itA = a.terms.begin();
        auto itB = b.terms.begin();

        while (itA != a.terms.end() && itB != b.terms.end()) {
            if (std::lexicographical_compare(itA->first.begin(), itA->first.end(),
                                             itB->first.begin(), itB->first.end()))
                return true;
            if (std::lexicographical_compare(itB->first.begin(), itB->first.end(),
                                             itA->first.begin(), itA->first.end()))
                return false;
            if (itA->second != itB->second)
                return itA->second < itB->second;
            ++itA;
            ++itB;
        }
        return a.terms.size() < b.terms.size();
    }
};

struct PairComparator {
    bool operator()(const std::pair<DifferentialForm, DifferentialForm>& a,
                    const std::pair<DifferentialForm, DifferentialForm>& b) const {
        Comparator cmp;
        if (cmp(a.first, b.first))
            return true;
        if (cmp(b.first, a.first))
            return false;
        return cmp(a.second, b.second);
    }
};