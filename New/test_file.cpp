#include "lib.h"

int main(int argc, char* argv[]) { // in the instruction, n can be 4 or 7
	// input the number n as reference for the N
	double n = atof(argv[1]); // convert char to double

	double N = pow(2,n);

	double gamma = 2; // implement w-cycles
	//double l = 8; //? as many multi-grid layers as possible

	double v1 = 1; // pre Gauss_Seidal iterations
	double v2 = 1; // post Gauss_Seidal iterations

	double h_f = 1 / N;

	// initialize U_0 with ZEROs
	// if N intervals and square, then (N+1) points per side
	std::vector<std::vector<double>> U(N+1, std::vector<double>(N+1, 0));
	//double res=func(1/N, 1/N);
	//std::cout << res << std::endl;
	std::vector<std::vector<double>> f(U.size(), std::vector<double>(U.size(), 0)); // Initialize RHS vector

	// The following is to initialize the f matrix and calculates the initial residual
	double r0 = 0;
	for (double i = 1; i < U.size()-1; i++) {
		for (double j = 1; j < U.size()-1; j++) {
			f[j][i] = func_RHS(j * h_f,(N-i) * h_f);
			if (abs(f[j][i]) > r0) {
				r0 = abs(f[j][i]);
			}
		}
	}

	double M = 20; // number of iterations for m(x_m, x_(m+1), x_(m+2)... etc)

	// initialize the U_out matrix
	//std::vector<std::vector<double>> u_out(N + 1, std::vector<double>(N + 1, 0));

	std::ofstream err_out;
	err_out.open("result_err_out.dat");

	for (int m = 1; m <= M; m++) {
		//std::vector<std::vector<double>> res = GS_smoother(U_0, f, v1);
		U = MG(U,/*U_0 is the initial guess to zero*/
			f, gamma, v1, v2);
		double max_residual=0;
		double res_ij=0;
		double residual_ratio;
		for (double i = 1; i < N; i++) {
			for (double j = 1; j < N; j++) {
				// res_ji is the residual for the point at (i,j)
				res_ij = f[j][i] + (U[j][i - 1] - 2 * U[j][i] + U[j][i + 1]) / (h_f * h_f) +
					(U[j - 1][i] - 2 * U[j][i] + U[j + 1][i]) / (h_f * h_f);

				if (abs(res_ij) > max_residual) {
					max_residual = abs(res_ij);
				}
			}
		}
		
		residual_ratio = max_residual / r0; // This is the residual ratio at m_th iteration

		// save the data
		err_out << m<< '\t' <<residual_ratio << std::endl;

		// show the data for checking
		std::cout << m << '\t' << residual_ratio << std::endl;
	}

	err_out.close();
	

	return 0;

}
