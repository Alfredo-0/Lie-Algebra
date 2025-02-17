//DifferentialForm.cpp

#include "Lie-Alg/DifferentialForm.h"

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
            for(auto i: indices){
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
    int sign = 1;
    
    std::array<int, DIMENSION> sortedIndices = indices;
    std::sort(sortedIndices.begin(), sortedIndices.begin()+degree);

    for (int i= 0; i< degree; ++i) {
        for (int j = i + 1; j < degree; ++j) {
            if (indices[i] > indices[j]) 
                sign = -sign;
        }
    }
    
    terms[sortedIndices] += (sign)*coeff;
}

DifferentialForm DifferentialForm::wedge(const DifferentialForm& other) const {
    int value = degree + other.degree;
    DifferentialForm result(value);
    bool duplicate = false;
    std::array<int, DIMENSION> combined;
    std::array<int, DIMENSION> count = {0};

    for (const auto& [indices1, coeff1] : this->terms) {
        for (auto& [indices2, coeff2] : other.terms) {

            for(int i = 0; i < DIMENSION; ++i){
                if (i < degree)
                    combined[i] = indices1[i];
                else if (i < degree + other.degree)
                    combined[i] = indices2[i-degree];
                else 
                    combined[i] = 0;
            }

            duplicate = false;
            count = {0};
            
            for (int& idx : combined) {
                if (++count[idx-1] == 2){
                    duplicate = true;
                    break; 
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
    std::array<int, DIMENSION> remaining = {};
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
            
            DifferentialForm remainingForm(degree - 1);

            double sign = ((j + 1) % 2 == 0) ? -1.0 : 1.0;
            
            remainingForm.addTerm(remaining, sign*coeff);
            
            result += remainingForm.wedge(algebra->dOf(indices[j])); 
            
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
            int i = 0;
            int j = 0;
            int k = 0;
            sign = 1.0;
            remaining = {0};
            while(i < degree && j < other.degree){
                if (indices[i] == indices2[j]){
                    aux = ((j+i) % 2 == 0) ? 1.0 : -1.0; 
                    sign *= aux;
                    ++j;
                    ++i;
                }
                else{
                    remaining[k]= indices[i];
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
                //result.print();
            }
        }
    }
    return result;
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
            for(auto i: indices){
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

LieAlgebra::LieAlgebra(std::vector<Pair> str) {
    std::cout << "Lie Algebra constructed.\n";
    
    for (int i = 0; i < DIMENSION; ++i) {
        DifferentialForm dEi(2);
        structureConstants[i] = dEi;
    }

    for(int it = 0; it < str.size(); it++){
        if(str[it].left != 0){
            structureConstants[it].addTerm({str[it].left/10, str[it].left%10}, str[it].right);
        }
    }
}
