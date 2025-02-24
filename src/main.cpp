//main.cpp
#include "Lie-Alg/DifferentialForm.h"
#include "Lie-Alg/PairUtils.h"
#include <iostream>

int main() {
    
    std::ifstream file("../../res/input.txt");
    if (!file) {
        std::cerr << "Error: Could not open input.txt for reading." << std::endl;
        return 1;
    }
    
    std::ofstream outfile("../../output.md");
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

        for(const auto &terms : lists.list2)
            for(const auto &p : terms)    
                omega.addTerm({p.left/10, p.left%10}, p.right);

        omega.inverse();

        outfile << "## Structure constants of the Lie Algebra:\n";
        outfile << "$(";

        for(const auto &terms : lists.list1){
            bool firstTerm = true;
            std::string index = " ";
            for(const auto &p : terms){
                if(!firstTerm)
                    outfile<<" + ";

                if(p.right == 0){
                    outfile << 0;
                    continue;
                }

                if(p.right<0)
                    outfile << '(' << p.right << ')';
                else if (p.right != 1)
                    outfile << p.right;

                outfile<<"e^{"<< p.left << " }";
                firstTerm = false;
            }
            outfile << ",\\ ";
        }
        outfile << ")$. \n\n";
        
        outfile << "### Symplectic form\n $\\omega=" << omega.toLaTeX()<<"$\n\n";
        
        if(omega.exteriorDerivative().checkZero()){
            std::cout<<"The symplectic form is closed!\n";
        }
        else{
            std::cout<<"The symplectic form is not closed!\n";
            omega.exteriorDerivative().print();
        }

        outfile<<"### Derivatives of $3-$forms\n";
        
        for(const auto& form : basis_3forms){
            DifferentialForm alpha(3);
            DifferentialForm dalpha(4);
            
            alpha.addTerm({form.i, form.j, form.k}, 1.0);
            dalpha = alpha.exteriorDerivative();

            if(!dalpha.checkZero()){
                outfile << "$d("<<alpha.toLaTeX() << ")  = " << dalpha.toLaTeX() << "$\n\n"; 
            }            
        }

        outfile<<"### Derivatives of $2-$forms\n";
        
        for(const auto& form : basis_2forms){
            DifferentialForm beta(2);
            DifferentialForm dbeta(3);
            
            beta.addTerm({form.i, form.j}, 1.0);
            dbeta = beta.exteriorDerivative();

            if(!dbeta.checkZero()){
                outfile << "$d("<<beta.toLaTeX() << ")  = " << dbeta.toLaTeX() << "$\n\n";
            }
            
        }   

        outfile<<"### $d \\Lambda d$ of $3-$forms\n";
        
        for(const auto& form : basis_3forms){
            DifferentialForm gamma(3);
            DifferentialForm dgamma(4);
            DifferentialForm ldgamma(2);
            DifferentialForm dldgamma(3);

            gamma.addTerm({form.i, form.j, form.k}, 1.0);
            
            dgamma = gamma.exteriorDerivative();
            ldgamma = dgamma.interiorProduct(omega.inverse());            
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
