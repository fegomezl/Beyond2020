#include "header.h"

double exact_w(const Vector &x)
{
  return 12*x(1)*x(1)-12*x(0)*x(0);
}
double exact_psi(const Vector &x)
{
  return pow(x(0),4)-pow(x(1),4);
}
