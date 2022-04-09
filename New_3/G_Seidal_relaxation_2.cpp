#include "lib_2.h"

void G_Seidal_smoother(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int v) {
	double h = 1 / double (u.size() - 1.);
	double h2 = h * h;
	std::vector<std::vector<double>> uold = u;

	for (int k = 0; k < v; k++) {
		uold = u;
		for (int i = 1; i < u.size()-1; i++) {
			for (int j = 1; j < u.size() - 1; j++) {
				u[i][j] = 0.25 * (h2 * f[i][j] + u[i - 1][j] + u[i][j - 1] + u[i][j + 1] + u[i + 1][j]);
			}
		}
	}


}