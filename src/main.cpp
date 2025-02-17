//main.cpp
#include <Lie-Alg/DifferentialForm.h>
#include "Lie-Alg/PairUtils.h"
#include <iostream>

int main() {
    PairLists lists = readPairLists("../res/input.txt");

    DifferentialForm::algebra = std::make_shared<LieAlgebra>(lists.list1);

    DifferentialForm omega(2);

    for(const auto &p : lists.list2)    
        omega.addTerm({p.left/10, p.left%10}, p.right);

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
    outfile << "### Symplectic form\n $\\omega=" << omega.toLaTeX()<<"$\n\n";
    
    outfile<<"### Derivatives\n";

    for(const auto& it : basis_3forms){
        DifferentialForm form(3);
        form.addTerm({it.i, it.j, it.k}, 1.0);
        DifferentialForm dform = form.exteriorDerivative();

        if(!dform.checkZero())
            outfile << "$d("<<form.toLaTeX() << ")  = " << dform.toLaTeX() << "$\n\n";             
    }
    
    for(const auto& it : basis_2forms){
        DifferentialForm form(2);
        form.addTerm({it.i, it.j}, 1.0);
        DifferentialForm dform = form.exteriorDerivative();

        if(!dform.checkZero())
            outfile << "$d("<<form.toLaTeX() << ")  = " << dform.toLaTeX() << "$\n\n";
        
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
        DifferentialForm d_form(4);
        DifferentialForm lambdad_form(2);
        DifferentialForm dlambdad_form(3);

        form.addTerm({it.i, it.j, it.k}, 1.0);
        
        d_form = form.exteriorDerivative();
        lambdad_form = d_form.interiorProduct(omega);
        dlambdad_form = lambdad_form.exteriorDerivative();

        if(!dlambdad_form.checkZero()){
            dlambdad_form.print();
            outfile << "$d \\Lambda d( " << form.toLaTeX() << ")  = " << dlambdad_form.toLaTeX() << "$\n\n";  
        }              
    }

    return 0;
}
