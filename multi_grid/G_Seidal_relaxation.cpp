#include "lib.h"

///////////////////////////////////////////////////////

double sq(double x) {
	return x * x;
}

double func(double x, double y) {
	//std::cout<< 8* M_PI*M_PI*sin(2*M_PI*x) * sin(2*M_PI*y)<<std::endl;	
	return 8* M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y); 
}

double U_analytical(double x, double y) { 
	return sin(2 * M_PI * x) * sin(2 * M_PI * y); 
}

std::vector<std::vector<double>> GS(std::vector<std::vector<double>> U_0, int iteration, double N) {
	double converge_threshold = 1e-10;

	std::vector<std::vector<double>> U(U_0);

	std::vector<std::vector<double>> U_old;
	double converged_max_error=0;

	double h = 1 / N; //mesh width because of N intervals

	for (int k = 0; k < iteration; k++) {

		U_old = U;
		double error_inf = 0;

		for (int i = 1; i < N; i++) {
			for (int j = 1; j < N; j++) {
				double func_result=func(i*h, j*h);	
				U[j][i] = 0.25 * (h*h * func_result + U[j][i - 1] + U[j][i + 1] + U[j - 1][i] + U[j + 1][i]);

				if (abs(U[j][i] - U_old[j][i]) > error_inf) {
					
					error_inf = abs(U[j][i] - U_old[j][i]);

				}
			}
		}

		std::cout <<"Error_inf in loop \t"<< error_inf << std::endl;

		if (error_inf < converge_threshold) {

			std::cout << k << std::endl;

			break;

		}


	}

	// measure the converged maximum error

	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {

			if (abs(U_analytical(i * h, j * h) - U[j][i])> converged_max_error) {
				double U_analytics=U_analytical(i*h, j*h);
				converged_max_error = abs(U_analytics - U[j][i]);

			}


		}
	}

	std::cout << "converged_max_error\t"<<converged_max_error << std::endl;

	return U;

	
}
// initialize u with zero
