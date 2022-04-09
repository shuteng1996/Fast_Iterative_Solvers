#include "lib.h"

int main(int argc, char* argv[]) {
	int N = 10;

	// initialize U_0 with ZEROs
	// if N intervals and square, then (N+1) points per side
	std::vector<std::vector<double>> U_0(N+1, std::vector<double>(N+1, 0));

	GS(U_0, 1000, N);

	return 0;

}