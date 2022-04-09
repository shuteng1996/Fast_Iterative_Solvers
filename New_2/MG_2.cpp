#include "lib_2.h"
std::vector<std::vector<double>> operator*(double val, std::vector<std::vector<double>> matrix) {
    std::vector<std::vector<double>> matrix_new(matrix.size(), std::vector<double>(matrix.size(), 0));
    for (double i = 0; i < matrix.size(); i++) {
        for (double j = 0; j < matrix.size(); j++) {
            matrix_new[j][i] = matrix[j][i] * val;
        }
    }

    return matrix_new;
}

std::vector<std::vector<double>> operator-(std::vector<std::vector<double>> matrix1, std::vector<std::vector<double>> matrix2) {
    if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
        std::cout << "warning, dimension not paired" << std::endl;
    }

    std::vector<std::vector<double>> result(matrix1.size(), std::vector<double>(matrix1.size(), 0));

    for (double i = 0; i < matrix1.size(); i++) {
        for (double j = 0; j < matrix1.size(); j++) {
            result[j][i] = matrix1[j][i] - matrix2[j][i];
        }
    }
    return result;

}



std::vector<std::vector<double>>Multigrid(std::vector<std::vector<double>> &u, std::vector<std::vector<double>>(&f), int gamma, int nu1, int nu2) {
    int N = u.size() - 1;
    int Nc = N / 2;
    double h = 1. / (double)N;
    double h2 = h * h;
    double H = 1. / (double)Nc;
    double H2 = H * H;

    smoother(u, f, nu1);
    std::vector<std::vector<double>>r = 0 * u;
    for (int i = 1; i < r.size() - 1; i++) {
        for (int j = 1; j < r[0].size() - 1; j++) {
            r[i][j] = f[i][j] + (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j] + u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / h2;
        }
    }
    std::vector<std::vector<double>>rc = RESTR(r);
    std::vector<std::vector<double>>ec = 0 * rc;
    std::vector<std::vector<double>>rcm = -1 * rc;
    if (ec.size() == 3) {
        //smoother(ec,rcm,1);
        ec[1][1] = 0.25 * H2 * rcm[1][1];
    }
    else {
        for (int i = 0; i < gamma; i++) {
            ec = Multigrid(ec, rcm, gamma, nu1, nu2);
        }
    }

    std::vector<std::vector<double>>el = PROLONG(ec);
    std::vector<std::vector<double>>temp = u - el;
    smoother(temp, f, nu2);
    return temp;
}
