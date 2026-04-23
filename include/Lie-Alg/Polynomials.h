#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <ginac/ginac.h>

extern GiNaC::symbol x;
extern GiNaC::ex expr;

// The letter I was changed for Z, since I is reserved to the imaginary unit.
extern GiNaC::symbol A, B, C, D, E, F, G, H;
extern GiNaC::symbol Z, J, K, L, M, N, O, P;
extern GiNaC::symbol Q, R, S, T;

extern std::map<std::string, GiNaC::symbol*> symbol_table;

extern GiNaC::ex A_hat, B_hat, C_hat, D_hat, E_hat, F_hat, G_hat, H_hat;
extern GiNaC::ex I_hat, J_hat, K_hat, L_hat, M_hat, N_hat, O_hat, P_hat;
extern GiNaC::ex Q_hat, R_hat, S_hat, T_hat;

extern std::map<std::string, GiNaC::ex*> expr_hat;

void initialize_symbols();
void extract_symbol_and_append(const GiNaC::ex& expr, GiNaC::lst& subs1, GiNaC::lst& subs3);

#endif // POLYNOMIALS_H
