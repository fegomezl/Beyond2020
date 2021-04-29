#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>


double initial(double r, double z, int m, int n);
double A(int m,int n);
double integrand(double r,double z,int n,int m);
double integrate(int m, int n);
double exact(double r, double z);


const double a=1.0;
const double b=1.0;

int main(void)
{


  return 0;
}

double initial(double r, double z, int m, int n)
{
return std::sin(M_PI*n*z/a)*std::cyl_bessel_j(0, boost::math::cyl_bessel_j_zero(0,m)*r/b);
}

double A(int m,int n)
{
  return 4*integrate(m,n)/(a*std::pow(b*std::cyl_bessel_j(1,boost::math::cyl_bessel_j_zero(0,m)),2));
}

double integrand(double r,double z,int n,int m)
{
  return initial(r,z,m,n)*std::sin(M_PI*n*z/a)*std::cyl_bessel_j(0, boost::math::cyl_bessel_j_zero(0,m)*r/b);
}

double integrate(int m, int n)
{ 
  auto f1 = [](double r, double z, int m, int n) { return integrand(r,z,m,n); };
  auto f = [&](double r, m) {
	     auto g = [&](double z, n) {
			return f1(r, z);
		      };
	     //return gauss_kronrod<double, 61>::integrate(g, 0, a, 5);
	     return boost::math::quadrature::trapezoidal(g, 0, a, 1e-6);
	   };
  double error;
  //double Q = gauss_kronrod<double, 15>::integrate(f, 0, b, 5, 1e-9, &error);
  double Q =  boost::math::quadrature::trapezoidal(f, 0, b, 1e-6);
  //std::cout << Q << ", error estimated at " << error <<std::endl;
  return Q;
}

double exact(double r, double z)
{
  double sum=0;
  
  for(int m=0; m<=0; m++)
    for(int n=0; n<=0; n++)
      sum+=A(m,n)*std::exp(1)*std::sin(1)*std::cyl_bessel_j(0,1);
  return sum;
}



