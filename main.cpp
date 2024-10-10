#include <iostream>
#include <iomanip>
#include <utility>
#include <cmath>
#include "brent.h"

/* Driver for routine brent */

#define TOL 1.0e-6

double func(double x) {
	return pow(x - 2.0, 2) - sqrt(2.0);
}

int main(void) {
	double xmin = 1.0, xmax = 3.0;

  std::cout << "Find minimum of f(x) = (x - 2)^2 - sqrt(2)" << std::endl;
	std::pair<double, double> res = brent(func, xmin, xmax, TOL);
	std::cout << std::setprecision(17) << "x_min, f_min = " << res.first << ", " << res.second << std::endl;

	return 0;
}