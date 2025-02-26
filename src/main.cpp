//main.cpp
#include "Lie-Alg/DifferentialForm.h"
#include "Lie-Alg/PairUtils.h"
#include <iostream>
#include <set>

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
    
    std::string line1, line2;
 
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

        outfile << "## Structure constants of the Lie Algebra:\n" << "$(";
        
        for(int it = 0; it < 6; ++it){
            if (omega.algebra->dOf(it).checkZero())
                outfile << 0;
            outfile << omega.algebra->dOf(it).toLaTeX();
            if (it != 5)
                outfile << ",\\ ";
        }
        outfile << ")$ \n\n";

        outfile << "### Symplectic form\n $\\omega=" << omega.toLaTeX() << "$\n\n"; 
        if(omega.exteriorDerivative().checkZero())
            std::cout<<"The symplectic form is closed!\n";
        else
            std::cout<<"The symplectic form is not closed!\n";

        std::multiset<DifferentialForm, Comparator> kernel;
        std::multiset<std::pair<DifferentialForm, DifferentialForm>, PairComparator> image;
        std::pair<DifferentialForm, DifferentialForm> pairForm;

        outfile<<"### Derivatives of $3-$forms\n";

        for(const auto& form : basis_3forms){
            DifferentialForm alpha({form.i, form.j, form.k}, 1.0);
            DifferentialForm dalpha = alpha.exteriorDerivative();

            if(!dalpha.checkZero()){
                pairForm = std::make_pair(dalpha, alpha);   
                image.insert(pairForm);
            }
            else{
                kernel.insert(alpha);
            }        
        }

        for(const auto& result: image)
            outfile << "$d("<<result.second.toLaTeX() << ")  = " << result.first.toLaTeX() << "$\n\n";

        outfile << "$Ker(d^3) \\supset \\{";
        for(const auto& result : kernel)
            outfile << result.toLaTeX() << ", \\ ";
        outfile << "\\}$ \n\n";

        image.clear();
        kernel.clear();

        outfile<<"### Derivatives of $2-$forms\n";

        for(const auto& form : basis_2forms){

            DifferentialForm beta({form.i, form.j}, 1.0);
            DifferentialForm dbeta = beta.exteriorDerivative();

            if(!dbeta.checkZero()){
                pairForm = std::make_pair(dbeta, beta);   
                image.insert(pairForm);
            }
            else{
                kernel.insert(beta);
            }
        }   

        for(const auto& result: image)
            outfile << "$d("<<result.second.toLaTeX() << ")  = " << result.first.toLaTeX() << "$\n\n";

        outfile << "$Ker(d^2) \\supset \\{";
        for(const auto& result : kernel)
            outfile << result.toLaTeX() << ", \\ ";
        outfile << "\\}$ \n\n";

        image.clear();
        kernel.clear();

        outfile<<"### $d \\Lambda d$ of $3-$forms\n";
        for(const auto& form : basis_3forms){
            DifferentialForm gamma({form.i, form.j, form.k}, 1.0);
            DifferentialForm dgamma = gamma.exteriorDerivative();
            DifferentialForm ldgamma = dgamma.interiorProduct(omega.inverse());
            DifferentialForm dldgamma = ldgamma.exteriorDerivative();

            if(!dldgamma.checkZero()){
                pairForm = std::make_pair(dldgamma, gamma);   
                image.insert(pairForm);
            } 
        }

        for(const auto& result: image)
            outfile << "$d("<<result.second.toLaTeX() << ")  = " << result.first.toLaTeX() << "$\n\n";
                 
        outfile << "\\pagebreak\n\n";
    }

    outfile.close();
    return 0;
}
