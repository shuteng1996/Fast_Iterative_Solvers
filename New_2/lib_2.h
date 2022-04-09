#ifndef LIB_2H_
#define LIB_2H_
#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>

// This defines the scalar-matrix multiplication
std::vector<std::vector<double>> operator*(double val, std::vector<std::vector<double>> matrix);

std::vector<std::vector<double>> operator-(std::vector<std::vector<double>> matrix1, std::vector<std::vector<double>> matrix2);

//std::vector<double> operator+(std::vector<double> vec1, std::vector<double> vec2);

//std::vector<double> operator-(std::vector<double> vec1, std::vector<double> vec2);

// This is the defined function in the instruction manual
double func_RHS(double x, double y);
double U_analytical(double x, double y);

std::vector<std::vector<double>> RESTR(std::vector<std::vector<double>> u);  // Restrict operator from fine grid to low grid
std::vector<std::vector<double>> PROLONG(std::vector<std::vector<double>> u); // Prolongation operator from coarse grid to find grid


void smoother(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int v);
std::vector<std::vector<double>>Multigrid(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>(&f), int gamma, int nu1, int nu2);

#endif // !LIB_2H_