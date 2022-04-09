#ifndef LIB_H_
#define LIB_H_
#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <vector>

double sq(double x);

// This is the defined function in the instruction manual
double func(double x, double y);

double U_analytical(double x, double y);



std::vector<std::vector<double>> GS(std::vector<std::vector<double>> U_0, int iteration, double N);


#endif // !LIB_H_
