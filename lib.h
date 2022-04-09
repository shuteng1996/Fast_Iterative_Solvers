#ifndef LIB_H_
#define LIB_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>

double sq(double x) {return x * x;}

// This is the defined function in the instruction manual
double func(double x, double y) { return 8 * sq(M_PI) * sin(2*M_PI*x)*sin(2*M_PI*y);}

double U_analytical(double x, double y) { return sin(2 * M_PI * x) * sin(2 * M_PI * y);}



std::vector<std::vector<double>> GS(std::vector<std::vector<double>> U_0, int iteration, int N);


#endif // !LIB_H_
