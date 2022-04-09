#include "lib_2.h"

double func_RHS(double x, double y) { return 8. * M_PI * M_PI * sin(2. * M_PI * x) * sin(2. * M_PI * y); }

double U_analytical(double x, double y) { return sin(2. * M_PI * x) * sin(2. * M_PI * y); }


int main(int argc, char* argv[]) {
    

    // prepare for the list of parameters needed
    int n = atoi(argv[1]);
    int M = atoi(argv[2]);

    int v1 = atoi(argv[3]);
    int v2 = atoi(argv[4]);
    int gamma = 2;

    int N = pow(2, n); // number of fine grid intervals
    int N_c = pow(2, n - 1);  // number of coarse grid intervals

    double h = 1. / (double) N; 


    std::vector<std::vector<double>>u(N + 1, std::vector<double>(N + 1, 0));
    std::vector<std::vector<double>>f(u.size(), std::vector<double>(u.size(), 0));
    std::vector<std::vector<double>>residual(u.size(), std::vector<double>(u.size(), 0));

    double r0 = 0; //initialize the residual

    for (int i = 1; i < f.size() - 1; i++) {
        for (int j = 0; j < f[0].size(); j++) {
            f[i][j] = func_RHS(i*h, j*h);
            if (abs(f[i][j] > r0)) { //f[i][j] is the residual here because of the zero U
                r0 = abs(f[i][j]); // update r0
            }

        }
    }


    std::ofstream err_output;
    err_output.open("err_2.dat");

    for (int m = 1; m <= M; m++) {
        u = MG(u, f, gamma, v1, v2); // looping the u, u_m, u_m+1, etc
        double max_residual = 0;
        double res_ij = 0;
        double residual_ratio;
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < N; j++) {
                // res_ji is the residual for the point at (i,j)
                res_ij = f[i][j] + (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j] + u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (h*h);
                
                if (abs(res_ij) > max_residual) {
                    max_residual = abs(res_ij); // update the max_residual for one total iteration
                }
            }
        }
        residual_ratio = max_residual / r0; // This is the residual ratio at m_th iteration
        //max_residual /= r0;
        err_output << m << '\t' << residual_ratio << std::endl;
        std::cout << "Iteration:\t " << m << "  Error:\t " << residual_ratio << std::endl;
    }
    err_output.close();





    return 0;
}
