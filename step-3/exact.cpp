#include <iostream>
#include <cmath>
#include <boost/math/quadrature/tanh_sinh.hpp>

double w(double x, double y);
double psi(double x, double y, double a, double b);
double integrand(double xp,double yp,double x,double y, double a, double b);
double Interate(double x, double y,double a,double b);

int main (void)
{
  double a=10.0, b=9.0;
  double dx=0.1, dy=0.1;
  double x=0, y=0;
  //exact solution for w!=0
  for(x=0; x<=a; x+=dx)
    for(y=0; y<=b; y+=dy)
      {
	std::cout << x <<"\t" << y <<"\t" << psi(x,y,a,b) << std::endl;
      }

  //exact solution for w=0

  return 0;
}
double w(double x, double y)
{
  return 1;
}
double psi(double x, double y, double a, double b)
{
  return Interate(x, y, a, b);
}
double integrand(double xp,double yp,double x,double y, double a, double b)
{		
  return std::log(std::hypot(x-xp,y-yp))*w(xp,yp);
}
double Interate(double x, double y, double a, double b) {

  boost::math::quadrature::tanh_sinh<double> integrator;
  auto inner_integral = [&integrator, x, y, a, b](double yp) {
        auto f = [yp, x, y, a, b](double xp) { return integrand(xp, yp, x, y, a, b); };
        return integrator.integrate(f, 0.0, a, 1e-8);
    };
    return integrator.integrate(inner_integral, 0.0, b, 1e-8);
}
