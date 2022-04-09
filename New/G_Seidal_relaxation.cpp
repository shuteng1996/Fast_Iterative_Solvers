#include "lib.h"

///////////////////////////////////////////////////////


//double sq(double x) {

	//return x * x;
//}

//double func_RHS(double x, double y) { 
	//return 8. * M_PI * M_PI * sin(2. * M_PI * x) * sin(2. * M_PI * y);
//}

//double U_analytical(double x, double y) { 
	//return sin(2 * M_PI * x) * sin(2 * M_PI * y); 
//}



// In this algorithm, U_0 and f are the FULL ones(including the boundary points)
std::vector<std::vector<double>> GS_smoother(std::vector<std::vector<double>> U_0, std::vector<std::vector<double>> f, int v) {
	double converge_threshold = 1e-10;
	double N = U_0.size() - 1;
	double h = 1 / N;

	std::vector<std::vector<double>> U(U_0);

	std::vector<std::vector<double>> U_old(U_0);

	//double converged_max_error=0;

	for (int k = 0; k < v; k++) {

		U_old = U;
		double error_inf = 0;
		double converged_max_error = 0;

		for (int i = 1; i < U_0.size()-1; i++) { // only work on the interior point values
			for (int j = 1; j < U_0.size()-1; j++) {

				//double func_result = func(i * h, j * h);

				U[j][i] = (1. / 4.) * (h * h * f[j][i] + U[j][i - 1] + U[j][i + 1] + U[j - 1][i] + U[j + 1][i]);

				//if (abs(U[j][i] - U_old[j][i]) > error_inf) {
					
					//error_inf = abs(U[j][i] - U_old[j][i]);

				//}

				//if (abs(U_analytical(j * h, (N-i) * h) - U[j][i]) > converged_max_error) {
					//double U_analytics = U_analytical(j * h, (N-i) * h);

					//converged_max_error = abs(U_analytics - U[j][i]);

				//}
			}
		}

		//std::cout <<"Error_inf in loop: \t"<< error_inf << std::endl;

		//std::cout << "converged_max_error\t" << converged_max_error << std::endl;
		//if (error_inf < converge_threshold) {

			//std::cout << k << std::endl;

			//break;

		//}


	}

	// measure the converged maximum error

	/*for (int i = 1; i < U_0.size()-1; i++) {
		for (int j = 1; j < U_0.size()-1; j++) {

			if (abs(U_analytical(j * h, (N-i) * h) - U[j][i])> converged_max_error) {
				double U_analytics = U_analytical(j * h, (N-i) * h);

				converged_max_error = abs(U_analytics - U[j][i]);

			}


		}
	}*/

	//std::cout << "converged_max_error\t"<<converged_max_error << std::endl;

	return U;

	
}
// initialize u with zero

//Test Gauss_Seidal

/*int main(int argc, char* argv[]) {
	double n = 4;
	double N = pow(2, n);
	double gamma = 2;
	double v1 = 30, v2 = 1;
	double hf = 1 / N;
	std::vector<std::vector<double>> U_0(N + 1, std::vector<double>(N + 1, 0));
	//double res=func(1/N, 1/N);
	//std::cout << res << std::endl;
	std::vector<std::vector<double>> f(U_0.size(), std::vector<double>(U_0.size(), 0)); // Initialize RHS vector
	double r0 = 0;
	for (double i = 1; i < U_0.size() - 1; i++) {
		for (double j = 1; j < U_0.size() - 1; j++) {
			f[j][i] = func_RHS(j * hf, (N - i) * hf);

		}
	}

	std::vector<std::vector<double>> result=GS_smoother(U_0, f, v1);
	return 0;

}*/

///
