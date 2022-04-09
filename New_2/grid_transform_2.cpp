#include "lib_2.h"


std::vector<std::vector<double>> RESTR(std::vector<std::vector<double>>u) {
    int Nc = (u.size() - 1) / 2;
    std::vector<std::vector<double>>uc(Nc + 1, std::vector<double>(Nc + 1, 0));
    for (int i = 1; i < Nc; i++) {
        int ii = 2 * i;
        for (int j = 1; j < Nc; j++) {
            int jj = 2 * j;
            uc[i][j] = 1. / 16. * (u[ii - 1][jj - 1] + 2 * u[ii][jj - 1] + u[ii + 1][jj - 1] + 2 * u[ii - 1][jj] + 4 * u[ii][jj] + 2 * u[ii + 1][jj] + u[ii - 1][jj + 1] + 2 * u[ii][jj + 1] + u[ii + 1][jj + 1]);
        }
    }
    return uc;
}

std::vector<std::vector<double>> PROLONG(std::vector<std::vector<double>> uc) {
    int Nc = uc.size() - 1;
    int N = 2 * Nc;
    std::vector<std::vector<double>>u(N + 1, std::vector<double>(N + 1, 0));
    for (int i = 1; i < Nc; i++) {
        int ii = 2 * i;
        for (int j = 1; j < Nc; j++) {
            int jj = 2 * j;
            u[ii][jj] += uc[i][j];
            u[ii + 1][jj] += 0.5 * uc[i][j];
            u[ii - 1][jj] += 0.5 * uc[i][j];
            u[ii][jj + 1] += 0.5 * uc[i][j];
            u[ii][jj - 1] += 0.5 * uc[i][j];
            u[ii + 1][jj + 1] += 0.25 * uc[i][j];
            u[ii + 1][jj - 1] += 0.25 * uc[i][j];
            u[ii - 1][jj + 1] += 0.25 * uc[i][j];
            u[ii - 1][jj - 1] += 0.25 * uc[i][j];
        }
    }
    return u;
}