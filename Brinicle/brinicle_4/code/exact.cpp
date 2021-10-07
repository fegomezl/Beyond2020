#include"header.h"

double initial_condition(double r, double z, int m, int n, double a, double b){
    return 0.;
}

double integrand(double r, double z, int m, int n, double a, double b){
    return 0.;
}

double integrate(int m, int n, double a, double b){
    auto inner_integral = [m,n,a,b](double z){
        auto f = [z,m,n,a,b](double r){
            return integrand(r,z,m,n,a,b); 
        };
        return 0.;
    };
    return 0.;
}

void Calc_Coe(double a, double b, std::vector<double> &Coeficients){
    for(int m = 1; m <= Mterms; m++)
        for(int n = 1; n <= Nterms; n++) 
            Coeficients[(m-1)+Mterms*(n-1)]=4*integrate(m,n,a,b)/(a*pow(1.,2));
}

double Aux(double r, double z, double t){
    int b = Rmax;
    int a = Zmax;
    double sum = 0;
    double Q;

    for(int m = 1; m <= Mterms; m++)
        for(int n = 1; n <= Nterms; n++){
            Q = Coeficients[(m-1)+Nterms*(n-1)];
            sum += 0.;
        }
    return sum;
}

double initial(const Vector &x, double t){
    return Aux(x(0),x(1),t) + 10.;
}
