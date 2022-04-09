#include "lib_2.h"


int main(int argc, char* argv[]) {
    double r0 = 0;
    int n = atoi(argv[1]);
    int M = atoi(argv[2]);
    int nu1 = atoi(argv[3]);
    int nu2 = atoi(argv[4]);
    int gamma = 2;

    //string pathname = "/home/felix/Dokumente/FIS/project2/MG_n" + to_string(n) + "_M" + to_string(M) + " _nu1" + to_string(nu1) + "_nu2" + to_string(nu2) + ".csv";

    int N = pow(2, n);
    int Nc = pow(2, n - 1);
    double h = 1. / (double)N;
    double h2 = h * h;

    std::vector<std::vector<double>>u(N + 1, std::vector<double>(N + 1, 0));
    std::vector<std::vector<double>>f(N + 1, std::vector<double>(N + 1, 0));
    std::vector<std::vector<double>>residual(N + 1, std::vector<double>(N + 1, 0));

    for (int i = 1; i < f.size() - 1; i++) {
        for (int j = 0; j < f[0].size(); j++) {
            f[i][j] = 8 * pow(M_PI, 2) * sin(2. * M_PI * i * h) * sin(2. * M_PI * j * h);
            if (abs(f[i][j] > r0)) {
                r0 = abs(f[i][j]);
            }

        }
    }
    std::ofstream mydata;
    mydata.open("err_2.dat");
    for (int m = 1; m <= M; m++) {
        u = Multigrid(u, f, gamma, nu1, nu2);
        double maxres = 0;
        double resij = 0;
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < N; j++) {
                resij = f[i][j] + (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j] + u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / h2;
                if (abs(resij) > maxres) {
                    maxres = abs(resij);
                }
            }
        }
        maxres /= r0;
        mydata << m << ";" << maxres << std::endl;
        std::cout << "Iteration: " << m << "  Error: " << maxres << std::endl;
    }
    mydata.close();





    return 0;
}
