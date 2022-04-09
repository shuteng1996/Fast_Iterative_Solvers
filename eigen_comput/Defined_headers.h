#ifndef DEFINED_HEADERS_H_
#define DEFINED_HEADERS_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "mmio.h"
#include <cmath>
#include <string>
#include <algorithm>
#include <utility>
#include <vector>
#include <fstream>
#include <chrono>

struct csr_matrix {
	std::vector<int> I;
	std::vector<int> J;
	std::vector<double> val;
	char sym = 'F'; // to tell if the matrix is symmetric or non-symmetric

	csr_matrix convert(std::string filename);
};

std::pair<std::vector<double>, std::vector<std::vector<double>>> Conjugate_gradient(csr_matrix A, std::vector<double> b,
	std::vector<double> x0, int iter, double tolerance_in_project);

std::vector<double> matrix_vector_unsymm(csr_matrix matrix, std::vector<double> x);
std::vector<double> matrix_vector_CSC(csr_matrix matrix, std::vector<double> x);
std::vector<double> matrix_vector_symm(csr_matrix matrix, std::vector<double> x);
std::vector<double> matrix_vector_prod(csr_matrix matrix_csr, std::vector<double> x);

double get_vector_norm(std::vector<double> vector);

double get_dot_product(std::vector<double> x, std::vector<double> y);

double get_norm_A(std::vector<double> x, csr_matrix A);

void get_Krylov_space(csr_matrix A, std::vector<std::vector<double>>(&V), std::vector<std::vector<double>>(&H), int j, int mycase);

std::pair<std::pair<std::vector<double>, double>,double> GMRES(csr_matrix A, std::vector<double> b, std::vector<double> x, int m, double tolerance_ratio_in_project, int mycase, std::vector<double> r0);

std::pair<std::vector<double>, double> rGMRES(csr_matrix A, std::vector<double> b, std::vector<double> x0, int m, double tolerance_ratio_in_project, int mycase);

template <typename T>
void swap_values(T& a, T& b);

std::vector<double> Preconditioning(csr_matrix A, std::vector<double> x, int case_num);

std::vector<double> operator/(std::vector<double> vec, double x);

#endif // !DEFINED_HEADERS_H_
