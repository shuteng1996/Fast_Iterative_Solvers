
#include "Defined_headers.h"





template <typename T> 
void swap_values(T &a, T &b) {
	T temp;
	temp = b;
	b = a;
	a = temp;
}


csr_matrix csr_matrix::convert(std::string filename) {
	int ret_code;
	MM_typecode matrix_code;
	// file pointer to read the file
	FILE* f = fopen(filename.c_str(), "r");
	int M, N, nz;  // M is the # of rows, N is the # of colums, nz is the # of elements
	int * I, * J;
	double* val;
	mm_read_banner(f, &matrix_code);
	mm_read_mtx_crd_size(f, &M, &N, &nz);

	// To tell if the matrix is symmetric or non-symmetric
	if (mm_is_symmetric(matrix_code)) {
		this->sym='T';
	}




	// in coordinate format, # of I, J and val are the same which are equal to nz
	I = (int*)malloc(nz * sizeof(int));
	J = (int*)malloc(nz * sizeof(int));
	val = (double*)malloc(nz * sizeof(double));
	
	// check for matrix dimensions
	std::cout << M <<' '<< N <<' '<< nz << '\n';
	std::cout << "Matrix Symmetric Properties: " << this->sym << '\n';

	for (int i = 0; i < nz; i++) {
		//%d means decimal interger, %lf means double type (Note: these values are not arranged in order)
		fscanf(f, "%d %d %lf\n", &I[i], &J[i], &val[i]);

		// To adapt to C matrix storage conventions, modify the matrix indexes
		I[i]--;
		J[i]--;
	}

	/* print out the values for testing */
	/*for (int j = 0; j < nz; j++) {
		printf("I:%d    J:%d      V:%f\n", I[j], J[j], val[j]);

	}
	*/
	/* rearrange the matrix in the oredr from low row number to high row number, for each row, from low column number
	to high column number*/

	for (int k = 0; k < nz-1; k++) {
		if (I[k] > I[k + 1]) {
			swap_values(I[k], I[k + 1]);
			swap_values(J[k], J[k + 1]);
			swap_values(val[k], val[k + 1]);

			// If this happens, then need to check the smaller value with the values before in a back order
			for (int i = k; i > 0; i--) {
				if (I[i] < I[i-1]) {
					swap_values(I[i], I[i - 1]);
					swap_values(J[i], J[i - 1]);
					swap_values(val[i], val[i - 1]);
				}
				else {
					break;
				}
			}
			//swap_values(I[k], I[k + 1]);
			//swap_values(J[k], J[k + 1]);
			//swap_values(val[k], val[k + 1]);
		}
	}

	/* for each row, check if the columns are arranged in increasing order */

	for (int k = 0; k < nz - 1; k++) {
		if (I[k] == I[k + 1] && J[k] > J[k + 1]) {
			swap_values(J[k], J[k + 1]);
			swap_values(val[k], val[k + 1]);
		}
	}

	/* print out the values for testing */
	std::cout << "After calibrating: " << std::endl;
	
	for (int j = 0; j < nz; j++) {
		printf("I:%d    J:%d      V:%f\n", I[j], J[j], val[j]);

	}
	


	/* Initialization with value zero */
	std::vector<int> num_elems_in_each_row(M,0);
	std::vector<int> newJ(nz ,0);
	std::vector<double> newVal(nz, 0);
	
	//std::vector<std::pair<double, double>> mypairs(nz);

	for (int i = 0; i < nz; i++) {
		num_elems_in_each_row[I[i]]++;

	}
	// num_elems_in_each_row before represents the true value

	for (int i = 1; i < M; i++) {
		num_elems_in_each_row[i] += num_elems_in_each_row[i - 1];
	}

	// use the iterator to do this

	num_elems_in_each_row.insert(num_elems_in_each_row.begin(), 0);

	//CSR format
	for (int k = 0; k < nz; k++) {
		newJ[k] = J[k];
		newVal[k] = val[k];
	}
	// test for CSR format
	/*
	std::cout << "After converting: " << std::endl;

	for (int j = 0; j < M+1; j++) {
		printf("I:%d\n", num_elems_in_each_row[j]);

	}
	*/
	this->I = num_elems_in_each_row;
	this->J = newJ;
	this->val = newVal;

	free(I);
	free(J);
	free(val);
	return *this;  // return the new struct

	// Up to now, the num_elems_in_each_row is index of I but in C notation
	// To get the real index, need to add 1 to each index

	/*std::vector<int> counter = num_elems_in_each_row;   */

	//for (int i = 0; i < nz; i++) {// nz here is the num of elems in each row
		//int u = I[i]; // the row of the ith data in C convention

		//assign the matrix in CSR format
		//newVal[counter[u]] = val[i];
		//newJ[counter[u]] = J[i];


		//mypairs[counter[u]] = std::make_pair(J[i], val[i]);
		/* mypairs show the value of J and val in C convention */
		//counter[u]++;
		/* counters give the corresponding J and val with respect to the given I */

		/* mypairs give the corresponding paired values for J and val as in given data, not necessarily in
		some order*/
	//}

	// This is to sort the matrix in ascending order by column
	//for (int i = 0; i < counter.size() - 1; i++) {
		//std::sort(mypairs.begin()+num_elems_in_each_row[i], mypairs.begin() + num_elems_in_each_row[i + 1], pair_compare);
	//}

	// assign the sorted values 
	
}

// matrix_vector (unsymmetrix version) (CSR matrix vector algorithm)
std::vector<double> matrix_vector_unsymm(csr_matrix matrix, std::vector<double> x) {
	// Here num_rows equals matrix.I.size()-1
	std::vector<double> y(matrix.I.size()-1,0);
	for (int i = 0; i < matrix.I.size() - 1; i++) {
		int i1 = matrix.I[i];
		int i2 = matrix.I[i + 1] - 1;
		for (int j = i1; j <= i2; j++) {
			y[i] = y[i] + matrix.val[j] * x[matrix.J[j]];
		}
	}

	return y;
}
// matrix vector (CSC matrix vector algorithm)
std::vector<double> matrix_vector_CSC(csr_matrix matrix, std::vector<double> x) {
	std::vector<double> y(matrix.I.size() - 1, 0);
	for (int i = 0; i < matrix.I.size() - 1; i++) {
		int i1 = matrix.I[i];
		int i2 = matrix.I[i + 1] - 1;
		for (int j = i1; j <= i2; j++) {
			y[matrix.J[j]] += matrix.val[j] * x[i];
		}
	}
	return y;
}

//matrix_vector product symmetric version
std::vector<double> matrix_vector_symm(csr_matrix matrix, std::vector<double> x) {
	std::vector<double> y1(matrix.I.size() - 1, 0);
	std::vector<double> y2(matrix.I.size() - 1, 0);
	std::vector<double> y_result(matrix.I.size() - 1, 0);


	y1 = matrix_vector_unsymm(matrix, x);
	y2 = matrix_vector_CSC(matrix, x);

	//print out the values for y2
	///for (auto i = y2.begin(); i != y2.end(); i++) {
		///std::cout << *i << std::endl;

	///}
	double diagonal_value;

	for (int i = 0; i < matrix.I.size() - 1; i++) {
		diagonal_value = 0;

		int i1 = matrix.I[i];
		int i2 = matrix.I[i + 1] - 1;
		for (int j = i1; j <= i2; j++) {
			if (matrix.J[j] == i) {
				diagonal_value = matrix.val[j];
				break;
			}
		}
		y_result[i] = y1[i] + y2[i] - x[i] * diagonal_value; // to counter the effect of considering the diagonal twice
	}
	return y_result;

}

std::vector<double> matrix_vector_prod(csr_matrix matrix_csr, std::vector<double> x) {
	std::vector<double> y_result;

	if (matrix_csr.sym == 'T') {
		y_result = matrix_vector_symm(matrix_csr, x);
	}
	else if (matrix_csr.sym == 'F') {
		y_result = matrix_vector_unsymm(matrix_csr, x);
	}

	///std::cout << "The final result of the matrix vector product is: " << '\n';
	///for (auto i = y_result.begin(); i != y_result.end(); i++) {
		///std::cout << *i << " ";
	///}
	///std::cout << '\n';

	return y_result;

}


