//main.cpp
#include "Lie-Alg/DifferentialForm.h"
#include "Lie-Alg/PairUtils.h"
#include <ginac/ginac.h>
#include <iostream>
#include <set>


int main() {

    GiNaC::symbol x("x", "\\lambda");

    std::ifstream file(RESOURCES_PATH);
    if (!file) {
        std::cerr << "Error: Could not open input.txt for reading." << std::endl << RESOURCES_PATH << std::endl;
        return 1;
    }
    
    std::ofstream outfile(DOCS_PATH);
    if (!outfile) {
        std::cerr << "Error: Could not open output.md for writing." << std::endl;
        return 1;
    }
    
    std::string enumeration, line1, line2;
 
    while(std::getline(file, enumeration, ' ')){
        if(enumeration.length() == 0)
            continue;
        
        
        std::getline(file, line1);
        std::getline(file, line2);
        
        PairLists lists = readPairLists(line1, line2, x);

        DifferentialForm::algebra = std::make_shared<LieAlgebra>(lists.list1);
        DifferentialForm omega(2);

        for(const auto &terms : lists.list2)
            for(const auto &p : terms)    
                omega.addTerm({p.left/10, p.left%10}, p.right);

        omega.inverse();

        outfile << "## " << enumeration << " Structure constants of the Lie Algebra:\n" << "$(";
        
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

        outfile<<"### $d \\Lambda d$ of $3-$forms\n";
        for(const auto& form : basis_3forms){
            DifferentialForm gamma({form[0], form[1], form[2]}, 1.0);
            DifferentialForm dgamma = gamma.exteriorDerivative();
            DifferentialForm ldgamma = dgamma.interiorProduct(omega.inverse());
            DifferentialForm dldgamma = ldgamma.exteriorDerivative();

            if(!dldgamma.checkZero()){
                pairForm = std::make_pair(dldgamma, gamma);   
                image.insert(pairForm);
            } 
        }

        for(const auto& result: image)
            outfile << "$" << result.second.getLetters() << ", \\ \\ d \\Lambda d("<<result.second.toLaTeX() << ")  = " << result.first.toLaTeX() << "\\ \\ " << result.first.getLetters() <<"$\n\n";


        image.clear();
        kernel.clear();
                 
        //outfile << "\\pagebreak\n\n";

    }

    outfile.close();
    return 0;
}