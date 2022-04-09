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



std::vector<std::vector<double>>MG(std::vector<std::vector<double>> &u, std::vector<std::vector<double>>(&f), int gamma, int v1, int v2) {
    int N = u.size() - 1;
    int N_c = N / 2;
    double h_f = 1. / (double)N;
    double h_c = 1. / (double)N_c;

    // perform G_Seidal smoothering operation
    G_Seidal_smoother(u, f, v1);
    std::vector<std::vector<double>>res(u.size(), std::vector<double>(u.size(),0));

    // calculate residuals
    for (int i = 1; i < res.size() - 1; i++) {
        for (int j = 1; j < res[0].size() - 1; j++) {
            res[i][j] = f[i][j] + (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j] + u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (h_f*h_f);
        }
    }

    // perform restriction operations
    std::vector<std::vector<double>>res_c = RESTR(res);



    std::vector<std::vector<double>>err_c = 0 * res_c;
    std::vector<std::vector<double>>minus_res_c = -1 * res_c;
    if (err_c.size() == 3) {
        // the coarest grid that can be achieved
        err_c[1][1] = 0.25 * h_c*h_c * minus_res_c[1][1];
    }
    else {
        for (int i = 0; i < gamma; i++) {// looping the MG algorithm

            err_c = MG(err_c, minus_res_c, gamma, v1, v2);
        }
    }

    std::vector<std::vector<double>>err = PROLONG(err_c);
    std::vector<std::vector<double>>corrected = u - err;
    G_Seidal_smoother(corrected, f, v2);
    return corrected;
}
