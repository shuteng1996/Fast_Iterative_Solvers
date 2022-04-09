// This is to implement Gram-Schmidt for the Krylov space to form an orthonormal basis

#include "Defined_headers.h"

inline double get_norm_A(std::vector<double> x, csr_matrix A) {
	return sqrt(get_dot_product(matrix_vector_prod(A, x), x));
}

double get_vector_norm(std::vector<double> vector) {
	double vector_norm_square=0;

	for (auto iter = vector.begin(); iter != vector.end(); iter++) {
		vector_norm_square = vector_norm_square + (*iter)*(*iter);
	}
	double vector_norm = sqrt(vector_norm_square);
	return vector_norm;

}

double get_dot_product(std::vector<double> x, std::vector<double> y) {


	double dot_product = 0;
	for (auto i = x.begin(),j = y.begin(); i != x.end(); i++, j++) {
		dot_product = dot_product + (*i) * (*j);
	}

	return dot_product;

}

// In the preconditioning section below, we assumed that the diagoals are non-zero

std::vector<double> Preconditioning(csr_matrix A, std::vector<double> x, int case_num) {
	// want case_num=0 to be Jacobin preconditioning, and case_num=1 to be G-S preconditioning

	// Note: x is the matrix to be preconditioned, y is the preconditione result
	int vector_size = x.size();
	
	std::vector<double> y(vector_size, 0);
	if (case_num == 0) {
		// No preconditioning processed
		y = x;
	}
	if (case_num == 1) {
		// Implement the Jacobi preconditioning
		
		for (int i = 0; i < vector_size; i++) {
			int i1 = A.I[i];
			int i2 = A.I[i + 1]-1;
			for (int j = i1; j <= i2; j++) {
				if (i == A.J[j]) {
					y[i] = x[i] / A.val[j];
					break;
				}
			}
		}
	}

	else if (case_num == 2) {
		//Implement G-S preconditioning
		for (int i = 0; i < A.I.size()-1; i++) {
			int i1 = A.I[i];
			int i2 = A.I[i + 1] - 1;
			for (int j = i1; j <= i2; j++) {
				if (i > A.J[j]) {
					// Then it is belongs to the lower triangular matrix including the diagonal
					x[i] -= A.val[j] * y[A.J[j]];
				}
				if (i == A.J[j]) {
					y[i] = x[i] / A.val[j];
				}

			}
		}
	}
	return y;
}

void get_Krylov_space(csr_matrix A, std::vector<std::vector<double>>(&V), std::vector<std::vector<double>>(&H), int j, int mycase/*mycase is for precondition*/) {
	// residual stands for the initial residual value 

	// Here starts the Gram-Schmidt process



	/////////////////////////
	// here just use the matrix_vector_prod which has been adjusted properly//
	/////////////////////////

	std::vector<double> w = matrix_vector_prod(A, V[j]);
	// Precondition Step:

	w = Preconditioning(A, w, mycase);

	for (int i = 0; i<=j; i++) {
		H[j][i] = get_dot_product(w, V[i]);
		// Over-writing step:
		for (int k = 0; k < w.size(); k++) {
			w[k] = w[k] - H[j][i] * V[i][k];
		}
			// Over_writing step done here.

	}
	H[j][j + 1] = get_vector_norm(w);
		
	// V[j+1] Over-writing step
	for (int k = 0; k < V[j + 1].size(); k++) {
		V[j + 1][k] = w[k] / H[j][j + 1];
	}

		// V[j+1] Over-writing step done
	


	// Here ends the Gram-Schmidt process

}

// redefine the minus operations for vector deduction

std::pair<std::pair<std::vector<double>,double>, double> GMRES(csr_matrix A, std::vector<double> b, std::vector<double> x, int m, double tolerance_ratio_in_project, int mycase/*mycase for precondition*/,
	std::vector<double> r0/*r0 is the preconditioned residual */) {
	std::ofstream res, krylov_ortho;

	res.open("error.dat", std::ofstream::app);
	krylov_ortho.open("krylov_orthogonal.dat", std::ofstream::app);  // output the orthogonality checking results

	// Initiaize the result x with the value of x0
	//std::vector<double> x(x);
	double error;

	// If the algorithm converge before m iterations, then record the stopping index
	int stopping_index=0;

	//A and x0 product, initialize with zero
	std::vector<double> r(b.size(),0);

	// initialize the Hessienberg matrix with m columns, each column consisting of m+1 values set to zero
	std::vector<std::vector<double>> H(m, std::vector<double>(m + 1, 0));

	//initialize the Karlov space with m columns, each column with the size of b.size() set to zero
	////To note here, V is a matrix of size m+1///////////////
	std::vector<std::vector<double >> V(m+1, std::vector<double>(b.size(), 0));

	// Initialize g values
	std::vector<double> g(m + 1, 0);

	std::vector<double> A_x_product = matrix_vector_prod(A, x);
	
	std::vector<double> C(m,0), S(m,0); // These represent the S and C in the equation, both having m values initialized to zero

	for (int i = 0; i < r.size(); i++) {
		r[i] = b[i] - A_x_product[i];
	}

	// Precondition Step:
	r=Preconditioning(A, r, mycase);

	double Beta = get_vector_norm(r);

	g[0] = Beta;// Initialize g[0] with Beta

	for (int i = 0; i < r.size(); i++) {
		V[0][i] = r[i] / (get_vector_norm(r));

	}
	
	//Invoke the get_Krylov space function
	for (int j = 0; j < m; j++) {
		get_Krylov_space(A, V, H, j, mycase);

		/////Test for Krylov space
		for (int i = 0; i < b.size(); i++) {
			std::cout << V[0][i] << "   ,    " << V[1][i] << std::endl;
		}
		std::cout << "Test for dot_product of V[0] and V[1]" << std::endl;
		std::cout << get_dot_product(V[0], V[1]) << std::endl;
		// Note: we ONLY need the Hessienberg matrix values from this procedure.

		for (int k = 1; k <= j; k++) {
			// H_j_k_1 is the copy of H[j][k-1]
			double H_j_k_1 = H[j][k - 1];
			// H_j_k is the copy of H[j][k]
			double H_j_k = H[j][k];

			H[j][k - 1] = C[k - 1] * H_j_k_1 + S[k - 1] * H_j_k;
			H[j][k] = -S[k - 1] * H_j_k_1 + C[k - 1] * H_j_k;
		}


		C[j] = H[j][j] / sqrt(H[j][j + 1] * H[j][j + 1] + H[j][j] * H[j][j]);
		S[j] = H[j][j + 1] / sqrt(H[j][j + 1] * H[j][j + 1] + H[j][j] * H[j][j]);

		H[j][j] = C[j] * H[j][j] + S[j] * H[j][j + 1];
		// After this step, can test that H[j][j+1] should be equal to zero.

		//The step for updating g[j+1] and g[j]
		g[j + 1] = -S[j] * g[j];
		g[j] = C[j] * g[j];

		error = abs(g[j + 1]); // For future output
		//Output the residual, orthogonality relations in a file

		res << abs(g[j + 1]) / get_vector_norm(r0) << std::endl;

		krylov_ortho << get_dot_product(V[0], V[j+1]) << std::endl;

		// record the j value and establish convergence criteria
		stopping_index=j;

		if ((abs(g[j + 1]) / get_vector_norm(r0)) < tolerance_ratio_in_project) {
			//std::cout << (stopping_index + 1) << std::endl;
			break;
		}
	}

	// end for
	
	//print out the stopping index
	//std::cout << stopping_index+1 << std::endl;



	// Initialize the resulting y vector

	//Initialize the y vector with size m and values zero.
	std::vector<double> y(m, 0);

	// Perform the back-substitution step
	for (int i = (stopping_index+1) - 1; i >= 0; i--) {
		for (int k = (stopping_index+1)-1 ; k > i; k--) {
			g[i] -= H[k][i]*y[k];
		}
		y[i] = g[i] / H[i][i];
	}

	// Back-substitution step done.

	//get the resulting x from x0
	for (int i = 0; i <= stopping_index; i++) {
		for (int j = 0; j < x.size(); j++) {
			x[j] =x[j] + y[i]*V[i][j];
		}
	}

	res.close();
	krylov_ortho.close();
	// Return the pair of final x_m and the residual
	return std::make_pair(std::make_pair(x, error),(stopping_index+1));

	


	
}

// Restarted GMRES algorithm
std::pair<std::vector<double>,double> rGMRES(csr_matrix A, std::vector<double> b, std::vector<double> x0, int m, double tolerance_ratio_in_project, int mycase) {
	
	// clear the original data
	std::ofstream res("error.dat", std::ofstream::trunc), krylov_ortho("krylov_orthogonal.dat", std::ofstream::trunc);

	
	//A and x0 product, initialize with zero
	double iteration_index=0;

	std::vector<double> r0(b.size(), 0);

	std::vector<double> A_xo_product = matrix_vector_prod(A, x0);



	for (int i = 0; i < r0.size(); i++) {
		r0[i] = b[i] - A_xo_product[i];
	}
	r0 = Preconditioning(A, r0, mycase);

	double rho = get_vector_norm(r0);

	std::vector<double> x;

	// initialize x with x0
	x = x0;

	while ((rho / get_vector_norm(r0)) >= tolerance_ratio_in_project) {

		auto result = GMRES(A, b, x, m, tolerance_ratio_in_project, mycase, r0);
		x = result.first.first;
		rho = result.first.second;

		iteration_index += result.second;


	}
	
	std::cout << iteration_index << std::endl;
	return std::make_pair(x, rho);

}