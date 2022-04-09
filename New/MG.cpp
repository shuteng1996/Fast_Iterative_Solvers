#include "lib.h"



std::vector<std::vector<double>> operator*(double val, std::vector<std::vector<double>> matrix) {
	std::vector<std::vector<double>> matrix_new(matrix.size(), std::vector<double>(matrix.size(),0));
	for (double i = 0; i < matrix.size(); i++) {
		for (double j = 0; j < matrix.size(); j++) {
			matrix_new[j][i] = matrix[j][i]*val;
		}
	}

	return matrix_new;
}

std::vector<std::vector<double>> operator-(std::vector<std::vector<double>> matrix1, std::vector<std::vector<double>> matrix2) {
	std::vector<std::vector<double>> result(matrix1.size(), std::vector<double>(matrix1.size(), 0));
	
	if (matrix1.size() !=matrix2.size() || matrix1[0].size()!=matrix2[0].size()) {
		std::cout<< "Warning, dimension not paired"<<std::endl;
	}

	for (double i = 0; i < matrix1.size(); i++) {
		for (double j = 0; j < matrix1.size(); j++) {
			result[j][i] = matrix1[j][i] - matrix2[j][i];
		}
	}
	return result;

}

/*std::vector<double> operator+(std::vector<double> vec1, std::vector<double> vec2) {
	std::vector<double> result(vec1.size(), 0);

	for (double i = 0; i < vec1.size(); i++) {
		result[i] = vec1[i] + vec2[i];
	}
	return result;

}*/

/*std::vector<double> operator-(std::vector<double> vec1, std::vector<double> vec2) {
	std::vector<double> result(vec1.size(), 0);

	for (double i = 0; i < vec1.size(); i++) {
		result[i] = vec1[i] - vec2[i];
	}
	return result;

}*/

double func_RHS(double x, double y) { return 8. * M_PI * M_PI * sin(2. * M_PI * x) * sin(2. * M_PI * y); }

double U_analytical(double x, double y) { return sin(2. * M_PI * x) * sin(2. * M_PI * y); }




std::vector<std::vector<double>> MG(std::vector<std::vector<double>> u, std::vector<std::vector<double>> f,
	double gamma, double v1, double v2) {

	double N = u.size() - 1; // fine grids
	double N_c = N / 2;  // coarse grids
	double h = 1 / N;  // fine grid width
	double h_c = 1 / N_c; // coarse grid width
	std::vector<std::vector<double>> u_smoothed = GS_smoother(u, f, v1); // v1 is the number pre-smoothing iterations

	// Initialize residual to zero
	std::vector<std::vector<double>> res(N + 1, std::vector<double>(N + 1, 0));

	// residual computation (Initialized here? )////
	for (double i = 1; i < N; i++) {
		for (double j = 1; j < N; j++) {

			res[j][i] = f[j][i] + (u_smoothed[j][i - 1] - 2 * u_smoothed[j][i] + u_smoothed[j][i + 1]) / (h * h) + (u_smoothed[j - 1][i] - 2 * u_smoothed[j][i] + u_smoothed[j + 1][i]) / (h * h);

		}
	}

	// finished residual computation

	// start restriction operation
	std::vector<std::vector<double>> res_c = RESTR(res);
	// end restriction operation

	// initialize error_coarse to zero
	std::vector<std::vector<double>> err_c(res_c.size(), std::vector<double>(res_c.size(), 0));


	if (err_c.size()==3) {
		err_c[1][1] =-1* (1. / 4.) * h_c * h_c * (res_c[1][1]);
	}
	else {
		
		for (double k = 0; k < gamma; k++) {
			//std::vector<std::vector<double>> minus_res=-1*res_c;
			err_c = MG( err_c, -1 * res_c, gamma, v1, v2);
		}
	}

	// started prolongation process
	std::vector<std::vector<double>> err = PROLONG(err_c);
	// finished prolongation process

	std::vector<std::vector<double>> u_out = GS_smoother(u_smoothed - err, f, v2);
	return u_out;
}
