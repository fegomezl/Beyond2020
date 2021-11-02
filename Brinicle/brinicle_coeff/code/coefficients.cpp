#include "header.h"

PhaseCoefficient::PhaseCoefficient(const GridFunction *Theta, const GridFunction *Phi, double v_l, double v_s, double inv_DT):
    Theta(Theta),
    Phi(Phi),
    v_l(v_l),
    v_s(v_s),
    inv_DT(inv_DT)
{}

double PhaseCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip){
    return 0.5*(1 + tanh(5*inv_DT*(Theta->GetValue(T, ip, 1)-T_fun(Phi->GetValue(T, ip, 1)))));
}
