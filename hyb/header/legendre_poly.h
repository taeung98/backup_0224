#include "main.h"

#ifndef L_POLY_H
#define L_POLY_H

void calculateLegendrePolynomial(int order, double x, double* result) {
  double p0 = 1.0;
  double p1 = x;
  *result = (order == 0) ? p0 : p1;

  for (int n = 2; n <= order; n++) {
    double pn = ((2 * n - 1) * x * p1 - (n - 1) * p0) / n;
    p0 = p1;
    p1 = pn;
    *result = pn;
  }
}

double P(int order, double value){ 

	double result;
	calculateLegendrePolynomial(order, value, &result);

  return result;
}

#endif
