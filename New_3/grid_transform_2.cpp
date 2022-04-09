#include "lib_2.h"


std::vector<std::vector<double>> RESTR(std::vector<std::vector<double>>u_f) {
    // Restrict operator from fine grid to coarse grid
    
    int N_c = (u_f.size() - 1) / 2;
    std::vector<std::vector<double>>u_c(N_c + 1, std::vector<double>(N_c + 1, 0));
    for (int i = 1; i < N_c; i++) { // coarse grid
        int ii = 2 * i;
        for (int j = 1; j < N_c; j++) {
            int jj = 2 * j;

            u_c[i][j] = 1. / 16. * (u_f[ii - 1][jj - 1] + 2 * u_f[ii][jj - 1] +
                u_f[ii + 1][jj - 1] + 2 * u_f[ii - 1][jj] + 4 * u_f[ii][jj] + 
                2 * u_f[ii + 1][jj] + u_f[ii - 1][jj + 1] + 2 * u_f[ii][jj + 1] + u_f[ii + 1][jj + 1]);
        }
    }
    return u_c; // coarse grid 
}

std::vector<std::vector<double>> PROLONG(std::vector<std::vector<double>> u_c) {
    int N_c = u_c.size() - 1;
    int N = 2 * N_c;

    std::vector<std::vector<double>>u_f(N + 1, std::vector<double>(N + 1, 0));
    
    for (int i = 1; i < N_c; i++) {
        int ii = 2 * i;
        for (int j = 1; j < N_c; j++) {
            int jj = 2 * j;
            u_f[ii][jj] += u_c[i][j];
            u_f[ii + 1][jj] += 0.5 * u_c[i][j];
            u_f[ii - 1][jj] += 0.5 * u_c[i][j];
            u_f[ii][jj + 1] += 0.5 * u_c[i][j];
            u_f[ii][jj - 1] += 0.5 * u_c[i][j];
            u_f[ii + 1][jj + 1] += 0.25 * u_c[i][j];
            u_f[ii + 1][jj - 1] += 0.25 * u_c[i][j];
            u_f[ii - 1][jj + 1] += 0.25 * u_c[i][j];
            u_f[ii - 1][jj - 1] += 0.25 * u_c[i][j];
        }
    }
    return u_f;
}