#include "lib.h"

std::vector<std::vector<double>> RESTR(std::vector<std::vector<double>> u_f) {
	// Restrict operator from fine grid to low grid
	double N_c = (1. / 2.) * (u_f.size() - 1); // coarse grid
	std::vector<std::vector<double>> u_c(N_c + 1, std::vector<double>(N_c + 1, 0)); // (N_c+1)*(N_c+1) matrix(coarse)
	double ii, jj;

	for (double i = 1; i < N_c; i++) {
		ii = 2. * i;
		for (double j = 1; j < N_c; j++) {
			jj = 2. * j;

			u_c[j][i] = (1. / 16.) * (u_f[jj - 1][ii - 1] + 2 * u_f[jj - 1][ii] + u_f[jj - 1][ii + 1]
				+ 2 * u_f[jj][ii - 1] + 4 * u_f[jj][ii] + 2 * u_f[jj][ii + 1]
				+ u_f[jj + 1][ii - 1] + 2 * u_f[jj + 1][ii] + u_f[jj + 1][ii + 1]);
		}
	}
	return u_c; // coarse grid
}

std::vector<std::vector<double>> PROLONG(std::vector<std::vector<double>> u_c) {
	double N_c = u_c.size() - 1; 
	double N = 2 * N_c;
	double ii, jj;

	std::vector<std::vector<double>> u_f(N+1, std::vector<double>(N+1, 0)); // initialize fine grid
	
	for (double i = 1; i < N_c; i++) {
		ii = 2 * i;
		for (double j = 1; j < N_c; j++) {
			jj = 2 * j;

			u_f[jj - 1][ii - 1] += (1. / 4.) * u_c[j][i];
			u_f[jj - 1][ii] += (1. / 2.) * u_c[j][i];
			u_f[jj - 1][ii + 1] += (1. / 4.) * u_c[j][i];
			u_f[jj][ii - 1] += (1. / 2.) * u_c[j][i];
			u_f[jj][ii] += u_c[j][i];
			u_f[jj][ii + 1] += (1. / 2.) * u_c[j][i];
			u_f[jj + 1][ii - 1] += (1. / 4.) * u_c[j][i];
			u_f[jj + 1][ii + 1] += (1. / 4.) * u_c[j][i];

		}
	}
	return u_f;


}