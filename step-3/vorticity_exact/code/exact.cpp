#include "header.h"

extern double exact_w(const Vector &x)
{
  return -8*pow(x(0),2)*x(1);
}
extern double exact_psi(const Vector &x)
{
  return x(1)*pow(x(0),4);
}
