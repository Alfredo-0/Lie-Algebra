//main.cpp
#include <Lie-Alg/DifferentialForm.h>
#include "Lie-Alg/PairUtils.h"
#include <iostream>

int main() {
    
    std::ifstream file("../res/input.txt");
    if (!file) {
        std::cerr << "Error: Could not open input.txt for reading." << std::endl;
        return 1;
    }
    
    std::ofstream outfile("../output.md");
    if (!outfile) {
        std::cerr << "Error: Could not open output.md for writing." << std::endl;
        return 1;
    }
    
    std::string line1;
    std::string line2;

    while(std::getline(file, line1)){
        
        if(line1.length() == 0)
            continue;
        std::getline(file, line2);
        
        PairLists lists = readPairLists(line1, line2);

        DifferentialForm::algebra = std::make_shared<LieAlgebra>(lists.list1);

        DifferentialForm omega(2);
        for(const auto &p : lists.list2)    
            omega.addTerm({p.left/10, p.left%10}, p.right);
        omega.print();
        
        outfile << "## Structure constants of the Lie Algebra:\n$";
        for (const auto &p : lists.list1){
            outfile << "(" << p.left<< ", " << p.right << ")";
        }
        outfile << "$\n\n";
        outfile << "### Symplectic form\n $\\omega=" << omega.toLaTeX()<<"$\n\n";
        
        outfile<<"### Derivatives\n";

        if(omega.exteriorDerivative().checkZero())
            std::cout<<"The symplectic form is closed!\n";
        else{
            std::cout<<"The symplectic form is not closed!\n";
            omega.exteriorDerivative().print();
        }

        for(const auto& it : basis_3forms){
            DifferentialForm alpha(3);
            alpha.addTerm({it.i, it.j, it.k}, 1.0);
            DifferentialForm dalpha = alpha.exteriorDerivative();

            if(!dalpha.checkZero()){
                outfile << "$d("<<alpha.toLaTeX() << ")  = " << dalpha.toLaTeX() << "$\n\n"; 
            }            
        }
        
        for(const auto& it : basis_2forms){
            DifferentialForm beta(2);
            beta.addTerm({it.i, it.j}, 1.0);
            DifferentialForm dbeta = beta.exteriorDerivative();

            if(!dbeta.checkZero()){
                outfile << "$d("<<beta.toLaTeX() << ")  = " << dbeta.toLaTeX() << "$\n\n";
            }
            
        }   
        
        for(const auto& it : basis_3forms){
            DifferentialForm gamma(3);
            DifferentialForm dgamma(4);
            DifferentialForm ldgamma(2);
            DifferentialForm dldgamma(3);

            gamma.addTerm({it.i, it.j, it.k}, 1.0);
            
            dgamma = gamma.exteriorDerivative();
            ldgamma = dgamma.interiorProduct(omega);            
            dldgamma = ldgamma.exteriorDerivative();

            if(!dldgamma.checkZero()){
                outfile << "$d \\Lambda d( " << gamma.toLaTeX() << ")  = " << dldgamma.toLaTeX() << "$\n\n";  
            }
        }
                 
        outfile << "\\pagebreak\n\n";
    }

    outfile.close();
    return 0;
}
