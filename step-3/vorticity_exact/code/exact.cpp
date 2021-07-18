#include "header.h"

double scaler=1e-4;

extern double exact_w(const Vector &x)
{
  return -8*pow(x(0),2)*x(1)*scaler;
}
extern double exact_psi(const Vector &x)
{
  return x(1)*pow(x(0),4)*scaler;
}
