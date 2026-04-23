#include "Lie-Alg/Polynomials.h"

// Definition of symbols
GiNaC::symbol x("x");
GiNaC::ex expr = x + x;

GiNaC::symbol A("A"), B("B"), C("C"), D("D"), E("E"), F("F"), G("G"), H("H");
GiNaC::symbol Z("Z", "I"), J("J"), K("K"), L("L"), M("M"), N("N"), O("O"), P("P");
GiNaC::symbol Q("Q"), R("R"), S("S"), T("T");

std::map<std::string, GiNaC::symbol*> symbol_table = {
        {"A", &A}, {"B", &B}, {"C", &C}, {"D", &D}, {"E", &E}, {"F", &F},
        {"G", &G}, {"H", &H}, {"Z", &Z}, {"J", &J}, {"K", &K}, {"L", &L},
        {"M", &M}, {"N", &N}, {"O", &O}, {"P", &P}, {"Q", &Q}, {"R", &R},
        {"S", &S}, {"T", &T}
};

GiNaC::ex A_hat, B_hat, C_hat, D_hat, E_hat, F_hat, G_hat, H_hat;
GiNaC::ex Z_hat, J_hat, K_hat, L_hat, M_hat, N_hat, O_hat, P_hat;
GiNaC::ex Q_hat, R_hat, S_hat, T_hat;

std::map<std::string, GiNaC::ex*> expr_hat = {
        {"A_hat", &A_hat}, {"B_hat", &B_hat}, {"C_hat", &C_hat}, {"D_hat", &D_hat}, {"E_hat", &E_hat}, {"F_hat", &F_hat},
        {"G_hat", &G_hat}, {"H_hat", &H_hat}, {"Z_hat", &Z_hat}, {"J_hat", &J_hat}, {"K_hat", &K_hat}, {"L_hat", &L_hat},
        {"M_hat", &M_hat}, {"N_hat", &N_hat}, {"O_hat", &O_hat}, {"P_hat", &P_hat}, {"Q_hat", &Q_hat}, {"R_hat", &R_hat},
        {"S_hat", &S_hat}, {"T_hat", &T_hat}
};

void extract_symbol_and_append(const GiNaC::ex& expr, GiNaC::lst& equations, GiNaC::lst& symbols, GiNaC::lst& coeff_values) {
//        const GiNaC::symbol* target_sym = nullptr;
//
//        for(size_t i = 0; i < expr.nops(); ++i) {
//                if (is_a<GiNaC::symbol>(expr.op(i))) {
//                        target_sym = &ex_to<GiNaC::symbol>(expr.op(i));
//                        equations.append(*target_sym == 0);
//                        symbols.append(*target_sym);
//                        if (i > 0 && is_a<GiNaC::numeric>(expr.op(i-1))){
//                                target_sym = &ex_to<GiNaC::symbol>(expr.op(i));
//                                coeff_values.append(*target_sym);
//                        }
//                        else{
//                                coeff_values.append(1);
//                        }
//                }
//                if (target_sym1 != nullptr) {
//                        std::cerr << "Warning: could not extract a symbol from coef1 = " << coef1 << std::endl;
//                        nope = false;
//                }
//                
//                
//        }
//
}

//// Try to extract a symbol from coef2
//const GiNaC::symbol *target_sym = nullptr;
//for (size_t i = 0; i < coef2.nops(); ++i) {
//    if (is_a<GiNaC::symbol>(coef2.op(i))) {
//        target_sym = &ex_to<GiNaC::symbol>(coef2.op(i));
//        break;
//    }
//}
//
//// If coef2 is itself a symbol
//if (is_a<GiNaC::symbol>(coef2)) {
//    target_sym = &ex_to<GiNaC::symbol>(coef2);
//}
//
//if (target_sym != nullptr) {
//    subs1.append(*target_sym == 0);
//    subs3.append(*target_sym);
//} else {
//    std::cerr << "Warning: could not extract a symbol from coef2 = " << coef2 << std::endl;
//    nope = false;
//}

////extract_symbol_and_append(coef2, sub1, sub3);
//
            //// Try to extract a symbol from coef1
            //const GiNaC::symbol *target_sym1 = nullptr;
            //for (size_t i = 0; i < coef1.nops(); ++i) {
            //    if (is_a<GiNaC::symbol>(coef1.op(i))) {
            //        target_sym1 = &ex_to<GiNaC::symbol>(coef1.op(i));
            //        break;
            //    }
            //}
//
            //// If coef1 is itself a symbol
            //if (is_a<GiNaC::symbol>(coef1)) {
            //    target_sym1 = &ex_to<GiNaC::symbol>(coef1);
            //}
//
            //if (target_sym1 != nullptr) {
            //    subs2.append(*target_sym1);
            //} else {
            //    std::cerr << "Warning: could not extract a symbol from coef1 = " << coef1 << std::endl;
            //    nope = false;
            //}


void initialize_symbols() {
    A_hat = A*(A*H - B*G - C*F - D*E + 2*Z*J + 2*K*L + 2*M*N - 2*O*P - 2*Q*R - 2*S*T)
            - 2*(B*(M*M - S*S) + C*(K*K - Q*Q) + E*(Z*Z - O*O) - B*C*E)
            - 4*(Z*K*M + O*K*S - Z*Q*S - O*Q*M);
        
    B_hat = B*(A*H - B*G + C*F + D*E + 2*Z*J + 2*K*L - 2*M*N - 2*O*P - 2*Q*R + 2*S*T)
            - 2*(-A*(N*N - T*T) + D*(K*K - Q*Q) + F*(Z*Z - O*O) + A*D*F)
            - 4*(Z*K*N + O*K*T - Z*Q*T - O*Q*N);

    C_hat = C*(A*H + B*G - C*F + D*E + 2*Z*J - 2*K*L + 2*M*N - 2*O*P + 2*Q*R - 2*S*T)
            - 2*(-A*(L*L - R*R) + D*(M*M - S*S) + G*(Z*Z - O*O) + A*D*G)
            - 4*(Z*L*M + O*L*S - Z*R*S - O*R*M);

    D_hat = D*(-A*H - B*G - C*F + D*E + 2*Z*J - 2*K*L - 2*M*N - 2*O*P + 2*Q*R + 2*S*T)
            - 2*(-B*(L*L - R*R) - C*(N*N - T*T) + H*(Z*Z - O*O) - B*C*H)
            - 4*(Z*L*N + O*L*T - Z*R*T - O*R*N);

    E_hat = E*(A*H + B*G + C*F - D*E - 2*Z*J + 2*K*L + 2*M*N + 2*O*P - 2*Q*R - 2*S*T)
            - 2*(-A*(J*J - P*P) + F*(M*M - S*S) + G*(K*K - Q*Q) + A*F*G)
            - 4*(J*K*M + P*K*S - J*Q*S - P*Q*M);

    F_hat = F*(-A*H - B*G + C*F - D*E - 2*Z*J + 2*K*L - 2*M*N + 2*O*P - 2*Q*R + 2*S*T)
            - 2*(-B*(J*J - P*P) + H*(K*K - Q*Q) - E*(N*N - T*T) - B*E*H)
            - 4*(J*K*N + P*K*T - J*Q*T - P*Q*N);

    G_hat = G*(-A*H + B*G - C*F - D*E - 2*Z*J - 2*K*L + 2*M*N + 2*O*P + 2*Q*R - 2*S*T)
            - 2*(-C*(J*J - P*P) - E*(L*L - R*R) + H*(M*M - S*S) - C*E*H)
            - 4*(J*L*M + P*L*S - J*R*S - P*R*M);

    H_hat = H*(-A*H + B*G + C*F + D*E - 2*Z*J - 2*K*L - 2*M*N + 2*O*P + 2*Q*R + 2*S*T)
            - 2*(-D*(J*J - P*P) - F*(L*L - R*R) - G*(N*N - T*T) + D*F*G)
            - 4*(J*L*N + P*L*T - J*R*T - P*R*N);

    Z_hat = Z*(A*H - B*G - C*F + D*E) - 2*J*(A*D - B*C)
            + 2*(A*(L*N - R*T) - B*(L*M - R*S) - C*(K*N - Q*T) + D*(K*M - Q*S))
            - 2*O*(Z*P - J*O + K*R - L*Q - M*T + N*S);

    J_hat = J*(-A*H + B*G + C*F - D*E) + 2*Z*(E*H - F*G)
            + 2*(E*(L*N - R*T) - F*(L*M - R*S) - G*(K*N - Q*T) + H*(K*M - Q*S))
            - 2*P*(Z*P - J*O + K*R - L*Q - M*T + N*S);

    K_hat = K*(A*H - B*G + C*F - D*E) - 2*L*(A*F - B*E)
            + 2*(A*(J*N + P*T) - B*(J*M + P*S) - E*(Z*N + O*T) + F*(Z*M + O*S))
            - 2*Q*(Z*P - J*O + K*R - L*Q + M*T - N*S);

    L_hat = L*(-A*H + B*G - C*F + D*E) + 2*K*(C*H - D*G)
            + 2*(C*(J*N + P*T) - D*(J*M + P*S) - G*(Z*N + O*T) + H*(Z*M + O*S))
            - 2*R*(Z*P - J*O + K*R - L*Q + M*T - N*S);

    M_hat = M*(A*H + B*G - C*F - D*E) - 2*N*(A*G - C*E)
            + 2*(A*(J*L - P*R) - C*(J*K - P*Q) - E*(Z*L - O*R) + G*(Z*K - O*Q))
            + 2*S*(Z*P - J*O - K*R + L*Q - M*T + N*S);

    N_hat = N*(-A*H - B*G + C*F + D*E) + 2*M*(B*H - D*F)
            + 2*(B*(J*L - P*R) - D*(J*K - P*Q) - F*(Z*L - O*R) + H*(Z*K - O*Q))
            + 2*T*(Z*P - J*O - K*R + L*Q - M*T + N*S);

    O_hat = O*(A*H - B*G - C*F + D*E) - 2*P*(A*D - B*C)
            - 2*(A*(L*T - R*N) - B*(L*S - R*M) - C*(K*T - Q*N) + D*(K*S - Q*M))
            - 2*Z*(Z*P - J*O + K*R - L*Q - M*T + N*S);

    P_hat = P*(-A*H + B*G + C*F - D*E) + 2*O*(E*H - F*G)
            - 2*(E*(L*T - R*N) - F*(L*S - R*M) - G*(K*T - Q*N) + H*(K*S - Q*M))
            - 2*J*(Z*P - J*O + K*R - L*Q - M*T + N*S);

    Q_hat = Q*(A*H - B*G + C*F - D*E) - 2*R*(A*F - B*E)
            + 2*(A*(J*T + P*N) - B*(J*S + P*M) - E*(Z*T + O*N) + F*(Z*S + O*M))
            - 2*K*(Z*P - J*O + K*R - L*Q + M*T - N*S);

    R_hat = R*(-A*H + B*G - C*F + D*E) + 2*Q*(C*H - D*G)
            + 2*(C*(J*T + P*N) - D*(J*S + P*M) - G*(Z*T + O*N) + H*(Z*S + O*M))
            - 2*L*(Z*P - J*O + K*R - L*Q + M*T - N*S);
            
    S_hat = S*(A*H + B*G - C*F - D*E) - 2*T*(A*G - C*E)
            + 2*(A*(J*R - P*L) - C*(J*Q - P*K) - E*(Z*R - O*L) + G*(Z*Q - O*K))
            + 2*M*(Z*P - J*O - K*R + L*Q - M*T + N*S);
            
    T_hat = T*(-A*H - B*G + C*F + D*E) + 2*S*(B*H - D*F)
            + 2*(B*(J*R - P*L) - D*(J*Q - P*K) - F*(Z*R - O*L) + H*(Z*Q - O*K))
            + 2*N*(Z*P - J*O - K*R + L*Q - M*T + N*S);
}

