//DifferentialForm.cpp
#include "Lie-Alg/DifferentialForm.h"

bool comp(int i, int j){
    if(i == 0 && j != 0)
        return false;
    if(i != 0 && j == 0)
        return true;
    return (i < j);
}

void DifferentialForm::print() const {
    std::string index = " ";
    bool firstTerm = true;

    for (auto& [indices, coeff]: terms){

        if(!firstTerm)
            std::cout<<" + ";

        if(coeff<0)
            std::cout << '(' << coeff << ')';
        else if (coeff > 0)
            std::cout << coeff;
        else
            continue;

        if(!indices.empty()){
            index = " ";
            for(auto const& i: indices){
                if(i == 0)
                    continue;
                index += std::to_string(i);
            }
            std::cout<<"e^{"<< index << " }";
        }
        firstTerm = false;
    }
    std::cout<<"\n"; 
}

bool DifferentialForm::checkZero() const{
    bool check = true;

    for (const auto& [indices, coeff] : this -> terms) {
        if (coeff != 0)
            return false;
    }
    return check;
}

void DifferentialForm::addTerm(const std::array<int, DIMENSION>& indices, double coeff) {
    if (!degree_assigned) {
        int computed_degree = 0;
        for (int idx : indices){
            if (idx != 0)
                ++computed_degree;
        }
        degree = computed_degree;
        degree_assigned = true;
    }
    
    double sign = 1.0;
    
    std::array<int, DIMENSION> sortedIndices = indices;
    std::sort(sortedIndices.begin(), sortedIndices.end(), comp);

    for (int i = 0; i < degree; ++i) {
        for (int j = i + 1; j < degree; ++j) {
            if (indices[i] > indices[j]) 
                sign = -sign;
        }
    }
    terms[sortedIndices] += (sign)*coeff;
    if(terms[sortedIndices] == 0)
        terms.erase(sortedIndices);
}

DifferentialForm DifferentialForm::wedge(const DifferentialForm& other) const {
    int value = degree + other.degree;
    bool duplicate = false;

    std::array<int, DIMENSION> combined = {0}, count = {0};
    DifferentialForm result(value);

    for (const auto& [indices1, coeff1] : this->terms) {
        for (auto& [indices2, coeff2] : other.terms) {
            duplicate = false;
            combined = {0};
            count = {0};

            for(int i = 0; i < value; ++i){
                if (i < degree){
                    combined[i] = indices1[i];
                    count[combined[i] - 1]++;
                }
                else if (i < degree + other.degree){
                    combined[i] = indices2[i-degree];
                    count[combined[i] - 1]++;
                }
                
                if (count[combined[i] - 1] > 1){
                    duplicate = true;
                    
                }
            }

            if (duplicate)
                continue;

            result.addTerm(combined, coeff1 * coeff2);
        }
    } 
    return result;
}

DifferentialForm& DifferentialForm::operator +=(const DifferentialForm& other){
    for(const auto& [indices, coeff] : other.terms)
        addTerm(indices, coeff);

    return *this;
}

DifferentialForm DifferentialForm::exteriorDerivative() const {
    DifferentialForm result(degree + 1);
    std::array<int, DIMENSION> remaining = {0};
    int aux = 0;

    for (const auto& [indices, coeff] : this->terms) {
        for (int j = 0; j < degree; ++j) {
            aux = 0;
            remaining = {0};

            for (int m = 0; m < degree; ++m) {
                if (m != j){
                    remaining[aux] = indices[m];
                    ++aux;
                }
            }
            
            double sign = ((j + 1) % 2 == 0) ? -1.0 : 1.0;
            
            DifferentialForm remainingForm(remaining, sign*coeff);
            
            result += remainingForm.wedge(algebra->dOf(indices[j] - 1)); 
        }
    }

    return result;
}

DifferentialForm DifferentialForm::interiorProduct(const DifferentialForm& other) const{
    DifferentialForm result;
    
    std::array<int, DIMENSION> remaining = {};
    double sign = 1.0;
    double aux = 0.0;

    for (const auto& [indices, coeff] : this->terms) { 
        for (const auto& [indices2, coeff2] : other.terms) {
            int i = 0, j = 0, k = 0;
            sign = 1.0;
            remaining = {0};
            while(i < degree && j < other.degree){
                if (indices[i] == indices2[j]){
                    aux = ((j + i) % 2 == 0) ? 1.0 : -1.0; 
                    sign *= aux;
                    ++j;
                    ++i;
                }
                else{
                    remaining[k] = indices[i];
                    ++k;
                    ++i;
                }
            }
            
            for (int id = i; id < degree; ++id){
                remaining[k] = indices[id];
                k++;
            }

            if (j == other.degree){
                result.addTerm(remaining, sign*coeff*coeff2);
            }
        }
    }
    return result;
}

DifferentialForm DifferentialForm::inverse() const {
    DifferentialForm inverse(2);

    for(auto const& [indices, coeff] : terms){
        inverse.addTerm({indices[0], indices[1]}, 1/coeff);
    }
    return inverse;
}

std::string DifferentialForm::toLaTeX() const {
    std::stringstream ss;
    std::string index = " ";
    bool firstTerm = true;

    for (auto& [indices, coeff]: terms){
        if(!firstTerm)
            ss<<" + ";

        if(coeff<0)
            ss << '(' << coeff << ')';
        else if (coeff != 1)
            ss << coeff;

        if(!indices.empty()){
            index = " ";
            for(auto const& i: indices){
                if(i == 0)
                    continue;
                index += std::to_string(i);
            }
            ss<<"e^{"<< index << " }";
        }
        firstTerm = false;
    }

    return ss.str();
}

LieAlgebra::LieAlgebra(std::vector<std::vector<Pair>> str) {
    std::cout << "Lie Algebra constructed.\n";
    
    for (int i = 0; i < 6; ++i) {
        DifferentialForm dEi(2);
        structureConstants[i] = dEi;
    }
    
    for(int it = 0; it < str.size(); it++){
        for(const auto& terms: str[it]){
            if(terms.left != 0){
                structureConstants[it].addTerm({terms.left/10, terms.left%10}, terms.right);
            }
        }
    }
}
