//main.cpp
#include <Lie-Alg/DifferentialForm.h>

#include <iostream>
#include "Lie-Alg/PairUtils.h"

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

int main() {

    PairLists lists = readPairLists("../res/input.txt");

    DifferentialForm::algebra = std::make_shared<LieAlgebra>(lists.list1);

    DifferentialForm omega(2);

    for(const auto &p : lists.list2)    
        omega.addTerm({p.left/10, p.left%10}, p.right);
    
    omega.print();

    std::ofstream outfile("../output.md");
    if (!outfile) {
        std::cerr << "Error: Could not open output.md for writing." << std::endl;
        return 1;
    }
    
    outfile << "## Structure constants of the Lie Algebra:\n$";
    for (const auto &p : lists.list1){
        outfile << "(" << p.left<< ", " << p.right << ")";
    }
    outfile << "$\n\n";
    outfile << "### Symplectic form\n $\\omega=\t";
    for (const auto &p : lists.list2){
        outfile <<p.right<< "\\cdot e^{" <<  p.left << "}+";
    }
    outfile<<"$\n\n";
    outfile<<"### Derivatives\n";

    for(const auto& it : basis_3forms){
        DifferentialForm form(3);

        form.addTerm({it.i, it.j, it.k}, 1.0);
        
        DifferentialForm dform = form.exteriorDerivative();

        if(!dform.checkZero()){
            std::string latexString = form.toLaTeX();
            std::string dlatexString = dform.toLaTeX();

            outfile << "$d("<<latexString<<")  = " << dlatexString << "$\n\n";  
        }              
    }
    
    for(const auto& it : basis_2forms){
        DifferentialForm form(2);

        form.addTerm({it.i, it.j}, 1.0);
        
        DifferentialForm dform = form.exteriorDerivative();

        if(!dform.checkZero()){
            std::string latexString = form.toLaTeX();
            std::string dlatexString = dform.toLaTeX();

            outfile << "$d("<<latexString<<")  = " << dlatexString << "$\n\n";
        }           
    }

    if(omega.exteriorDerivative().checkZero())
        std::cout<<"The symplectic form is closed!\n";
    else{
        std::cout<<"The symplectic form is not closed!\n";
        omega.print();
        omega.exteriorDerivative().print();
    }
    
    for(const auto& it : basis_3forms){
        DifferentialForm form(3);

        form.addTerm({it.i, it.j, it.k}, 1.0);
        
        DifferentialForm d_lambda_dform = form.exteriorDerivative().interiorProduct(omega).exteriorDerivative();

        if(!d_lambda_dform.checkZero()){
            d_lambda_dform.print();
            std::string latexString = form.toLaTeX();
            std::string dlatexString = d_lambda_dform.toLaTeX();

            outfile << "$d \\Lambda d( "<<latexString<<")  = " << dlatexString << "$\n\n";  
        }              
    }

    DifferentialForm alpha(3);
    alpha.addTerm({2, 4, 5}, 1.0);
    alpha.print();
    
    DifferentialForm dalpha(4);
    dalpha = alpha.exteriorDerivative();
    dalpha.print();
    
    DifferentialForm lambdaalpha(2);
    lambdaalpha = dalpha.interiorProduct(omega);
    lambdaalpha.print(); 
    
    lambdaalpha.exteriorDerivative().print(); 

    return 0;
}
