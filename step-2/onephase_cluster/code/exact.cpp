#include"header.h"

double initial_condition(double r, double z, int m, int n, double a, double b)
{
  return r*sin(n*M_PI*z/a)*boost::math::cyl_bessel_j(0,boost::math::cyl_bessel_j_zero(0.0,m)*r/b);
}

void Calc_Coe(double a, double b, std::vector<double> & Coeficients)
{
  for(int m=1; m<=Mterms; m++)
    for(int n=1; n<=Nterms; n++) 
        Coeficients[(m-1)+Mterms*(n-1)]=4*integrate(m,n,a,b)/(a*pow(b*boost::math::cyl_bessel_j(1.0, boost::math::cyl_bessel_j_zero(0.0,m)),2));
}


double integrand(double r,double z,int m,int n, double a, double b)
{
  return initial_condition(r,z,1,1,a,b)*sin(M_PI*n*z/a)*boost::math::cyl_bessel_j(0.0,boost::math::cyl_bessel_j_zero(0.0,m)*r/b)*r;
}


double integrate(int m, int n,double a,double b)
{

    auto inner_integral = [m,n,a,b](double z)
        {
          auto f = [z,m,n,a,b](double r) { return integrand(r,z,m,n,a,b);} ;
          return boost::math::quadrature::trapezoidal(f,0.0,a,1e-8);};
    return boost::math::quadrature::trapezoidal(inner_integral,0.0,b,1e-8);
}


double initial(const Vector &x, double t)
{
  return Aux(x(0),x(1),t);
}

void print_exact()
{

	ofstream myfile;
    myfile.open ("results/initial_exact.txt");
    double r, z;
    for(r=0; r<=Rmax; r+=0.1){
        for(z=0; z<=Zmax; z+=0.1 ){
        	myfile << r << "\t" << z << "\t" << Aux(r, z, 0) << endl;}}
     myfile.close();
}
void print_initial()
{
	ofstream myfile;
    myfile.open ("results/initial.txt");
    double r, z;
    for(r=0; r<=Rmax; r+=0.1){
        for(z=0; z<=Zmax; z+=0.1 ){
        	myfile << r << "\t" << z << "\t" << initial_condition( r, z, 1, 1, Zmax, Rmax) << endl;}}
     myfile.close();
}

double Aux( double r, double z, double t)
{
  int b=Rmax;
  int a=Zmax;
  double alpha=alpha_l;
  double sum=0;
  double Q;

  for(int m=1; m<=Mterms; m++)
    for(int n=1; n<=Nterms; n++)
    {
      Q=Coeficients[(m-1)+Nterms*(n-1)];
      sum +=Q*exp(-(pow((n*M_PI/a),2)+pow((boost::math::cyl_bessel_j_zero(0.0,m)/b),2))*alpha*t)*sin(n*M_PI*z/a)*boost::math::cyl_bessel_j(0.0,boost::math::cyl_bessel_j_zero(0.0,m)*r/b);
    }

  return sum;
}

