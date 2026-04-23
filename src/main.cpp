//main.cpp
#include "Lie-Alg/DifferentialForm.h"
#include "Lie-Alg/Polynomials.h"
#include "Lie-Alg/PairUtils.h"
#include <ginac/ginac.h>
#include <iostream>
#include <set>


int main() {
    initialize_symbols();

    GiNaC::symtab table;
    for (const auto& [name, sym_ptr] : symbol_table) {
        table[name] = *sym_ptr;
    }

    GiNaC::parser reader(table);

    GiNaC::symbol x("x");

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

        // new code here
        outfile<<"### Derivatives of $3-$forms\n";

        for(const auto& terms : primitive_basis_3forms){
            DifferentialForm gamma;
            if (terms[0] == terms[1])
                gamma.addTerm({basis_3forms[terms[0]][0],basis_3forms[terms[0]][1],basis_3forms[terms[0]][2]}, 1.0);
            else if (terms[0] < terms[1]){
                gamma.addTerm({basis_3forms[terms[0]][0],basis_3forms[terms[0]][1],basis_3forms[terms[0]][2]}, 1.0);
                gamma.addTerm({basis_3forms[terms[1]][0],basis_3forms[terms[1]][1],basis_3forms[terms[1]][2]}, -1.0);
            }
            else{
                gamma.addTerm({basis_3forms[terms[0]][0],basis_3forms[terms[0]][1],basis_3forms[terms[0]][2]}, 1.0);
                gamma.addTerm({basis_3forms[terms[1]][0],basis_3forms[terms[1]][1],basis_3forms[terms[1]][2]}, 1.0);
            }

            DifferentialForm dgamma = gamma.exteriorDerivative();
            
            if(!dgamma.checkZero()){
                pairForm = std::make_pair(dgamma, gamma);   
                image.insert(pairForm);
            }         
        }

        GiNaC::lst closed_conditions;

        for(const auto& result: image){
            GiNaC::ex coef1 = reader(result.second.getLetters());
            outfile << "$" << coef1 <<", \\ \\ d("<<result.second.toLaTeX() << ")  = " << result.first.toLaTeX() << "$\n\n";
            closed_conditions.append(coef1 == 0); 
        }

        image.clear();
        kernel.clear();
        
        // new code ends here

        outfile<<"### $d \\Lambda d$ of $3-$forms.\n\n";

        for(const auto& terms : primitive_basis_3forms){
            DifferentialForm gamma;
            if (terms[0] == terms[1])
                gamma.addTerm({basis_3forms[terms[0]][0],basis_3forms[terms[0]][1],basis_3forms[terms[0]][2]}, 1.0);
            else if (terms[0] < terms[1]){
                gamma.addTerm({basis_3forms[terms[0]][0],basis_3forms[terms[0]][1],basis_3forms[terms[0]][2]}, 1.0);
                gamma.addTerm({basis_3forms[terms[1]][0],basis_3forms[terms[1]][1],basis_3forms[terms[1]][2]}, -1.0);
            }
            else{
                gamma.addTerm({basis_3forms[terms[0]][0],basis_3forms[terms[0]][1],basis_3forms[terms[0]][2]}, 1.0);
                gamma.addTerm({basis_3forms[terms[1]][0],basis_3forms[terms[1]][1],basis_3forms[terms[1]][2]}, 1.0);
            }

            DifferentialForm dgamma = gamma.exteriorDerivative();
            DifferentialForm ldgamma = dgamma.interiorProduct(omega.inverse());
            DifferentialForm dldgamma = ldgamma.exteriorDerivative();

            if(!dldgamma.checkZero()){
                pairForm = std::make_pair(dldgamma, gamma);   
                image.insert(pairForm);
            } 

        }
        
        GiNaC::lst vars;
        GiNaC::lst subs1;
        GiNaC::lst subs2;
        GiNaC::lst equations;
        
        std::set<GiNaC::symbol, GiNaC::ex_is_less> symbol_set;
        
        for(const auto& result: image){
            GiNaC::ex coef1 = reader(result.second.getLetters());
            GiNaC::ex coef2 = reader(result.first.getLetters());
            
            std::string name = ex_to<GiNaC::symbol>(coef1).get_name() + "_hat";  
            subs2.append(coef1 == *expr_hat[name]);
            
            equations.append(coef2 == coef1); 

            if ( coef2.nops() == 0 ){
                symbol_set.insert(ex_to<GiNaC::symbol>(coef2));
            }
            else{
                for (size_t i = 0; i < coef2.nops(); ++i) {
                    if(coef2.op(i).nops() > 0){
                        symbol_set.insert(ex_to<GiNaC::symbol>(coef2.op(i).op(0)));
                    }
                    else if (is_a<GiNaC::symbol>(coef2.op(i))) {
                        symbol_set.insert(ex_to<GiNaC::symbol>(coef2.op(i)));
                    }
                }
            }
            outfile << "$" << coef1 << ", \\ \\ d \\Lambda d("<<result.second.toLaTeX() << ")  = " << result.first.toLaTeX() << "\\ \\ " << coef2 <<"$\n\n";
        }
        
        for (const GiNaC::symbol &s : symbol_set) {
            vars.append(s);
            subs1.append(s == 0);
        }
        
        GiNaC::lst solutions = ex_to<GiNaC::lst>(GiNaC::lsolve(equations, vars));

        std::cout<<equations<<" - "<<solutions<<std::endl;

        outfile << "### ODE system \n";

        outfile << GiNaC::latex;

        int counter = 0;
        for (const auto sol : solutions){
            outfile << "$\\partial_t "<< sol.lhs() << " = " << collect((sol.rhs().subs(subs2)).expand() - (sol.rhs().subs(subs2).subs(subs1)).expand(), vars) << "$\n" <<std::endl;
            counter++;    
        }

        outfile << "### ODE system with closed 3-forms\n";

        counter = 0;
        for (const auto sol : solutions){
            outfile << "$\\partial_t "<< sol.lhs() << " = " << collect((sol.rhs().subs(subs2)).expand() - (sol.rhs().subs(subs2).subs(subs1)).expand(), vars).subs(closed_conditions) << "$\n" <<std::endl;
            counter++;    
        }

        

        image.clear();
        kernel.clear();
                 
        outfile << "\\pagebreak\n\n";
    }

    outfile.close();
    return 0;
}