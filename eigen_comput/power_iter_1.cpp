#include "Defined_headers.h"

inline std::vector<double> operator/(std::vector<double> vec, double x) {
	   std::vector<double> result(vec.size(),0);
	   for (int i = 0; i != vec.size()-1; i++) {
		   result[i] = vec[i] / x;
	   }
	   return result;
}

inline std::vector<double> operator*(double scalar, std::vector<double> vec) {
	std::vector<double> result(vec.size(), 0);
	for (int i = 0; i != vec.size() - 1; i++) {
		result[i] = vec[i] * scalar;
	}
	return result;
}

inline std::vector<double> operator-(std::vector<double> vec1, std::vector<double> vec2) {
	std::vector<double> result(vec1.size(), 0);
	for (int i = 0; i != vec1.size() - 1; i++) {
		result[i] = vec1[i] - vec2[i];
	}
	return result;
}


double power_iter(csr_matrix A, double threshold=1e-8) {
	// return the size of the matrix
	double N = A.I.size() - 1;
	//initialize q0 
	std::vector<double> q(N, 1/sqrt(N));
	//declare lambda
	double lambda_old = 0.;
	double lambda = 0.;
	int iteration = 0.;

	// initialize the error
	double error = 1.;
	std::vector<double> z = matrix_vector_prod(A, q);
	std::ofstream output("result.dat");

	// start timing
	auto time_start = std::chrono::high_resolution_clock::now();





	///////////////////////////////////////////
	// looping
	while (error >= threshold) {
		lambda_old = lambda;

		q = z / get_vector_norm(z);

		z = matrix_vector_prod(A, q);// z is calculated here for convenience to get lambda
		lambda = get_dot_product(q, z);
				

		output<< lambda<< std::endl;
		

		error = abs(lambda - lambda_old);
		iteration++;
 		//output<<error<<std::endl;
		std::cout << "iteration index" << '\t' << iteration << '\t' << error << std::endl;

	}


	auto time_end = std::chrono::high_resolution_clock::now();

	output.close();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_start);


	//////////////////////////////////////////////
	std::cout << "Elapsed time(s): " << '\t' << elapsed.count() * 1e-9 << std::endl;
	return lambda;
}

//Step 2: For the Lancoz algorithm
double power_iter_Lanczos(std::vector<double> (&alphas), std::vector<double> (&betas), double threshold, int m) {
	std::vector<double> q(m, 1. / sqrt((double)m));

	double theta_old = 0.; // theta is the eigenvalue of the Lanczos matrix
	double theta = 0.;
	int iteration = 0;
	double err = 1.0;

	std::vector<double> z(q);
	

	std::ofstream output1("lanczos_result.dat");



	z[0] = alphas[0] * q[0] + betas[0] * q[1];

	for (int i = 1; i < m - 1; i++) {
		z[i] = betas[i - 1] * q[i - 1] + alphas[i] * q[i] + betas[i + 1] * q[i + 1];
	}
	z[m - 1] = betas[m - 2] * q[m - 2] + alphas[m - 1] * q[m - 1];

	while (err >= threshold) {
		theta_old = theta;

		q = z / get_vector_norm(z);

		// calculate the next z
		z[0] = alphas[0] * q[0] + betas[0] * q[1];
		for (int i = 1; i < m-1; i++) {
			z[i] = betas[i - 1] * q[i - 1] + alphas[i] * q[i] + betas[i + 1] * q[i + 1];
		}
		z[m - 1] = betas[m - 2] * q[m - 2] + alphas[m - 1] * q[m - 1]; /////////////////
		// At this step, next z is calculated
		theta = get_dot_product(q, z);
		err = abs(theta_old - theta);
		
		
		output1<< theta<< std::endl;


		iteration++;

		std::cout << "Lanczos matrix iteration index" << '\t' << iteration << '\t' << err << std::endl;
	}

	output1.close();
	return theta;
}



/* Lanczos algorithm is used for symmetric matrices*/
double lanczos(csr_matrix A, int m, double threshold) {
	   // returns the matrix dimension
	   double N = A.I.size() - 1;

	   std::vector<double> w;

	   //initialize the first v
	   std::vector<double> v(N, 1./ sqrt((double)N));

	   std::vector<double> vold(N, 0);
	   // initialize the diagonal and off-diagonal values
	   std::vector<double> alpha_vector, beta_vector;

	   // initial alpha, beta
	   double alpha = 0, beta = 0;

	   // start timing
	   auto time_start=std::chrono::high_resolution_clock::now();



	   /////////////////////////////////////
	   for (int i = 0; i < m; i++) {
		   w = matrix_vector_unsymm(A, v);
		   w=w-beta*vold;

		   alpha = get_dot_product(w, v);
		   w = w - alpha * v;

		   beta = get_vector_norm(w);

		   vold = v;
		   v = w / beta;

		   alpha_vector.push_back(alpha);
		
		   beta_vector.push_back(beta); // Notice: should take out the last one
	
	   }

	   // Finished converting to the m*m matrix
	   double eigen_val=power_iter_Lanczos(alpha_vector, beta_vector, threshold, m);

	   // record the time
	   auto time_end = std::chrono::high_resolution_clock::now();

	   auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_start);


	   //////////////////////////////////////////////
	   std::cout << "Elapsed time(s): "<<'\t'<< elapsed.count() * 1e-9 << std::endl;

	   return eigen_val;
	   // timing

}

int main(int argc, char* argv[]) {

	int m; // the size of the Lanczos matrix
	double threshold; // threshold value
	double eigenvalue;

	csr_matrix A;

	int implem_choice;
	implem_choice = atoi(argv[1]);

	if (implem_choice == 0) {// This means applying the power iteration directly to nos6.mtx
		// argv[2] is the matrix name
		A = A.convert(argv[2]);
		eigenvalue=power_iter(A);
		std::cout << " nos6.mtx: eigenvalue: " << '\t' << eigenvalue << std::endl;
	}
	else if (implem_choice == 1) { // This means applying the power iteration directly to s3rmt3m3.mtx
		A = A.convert(argv[2]);
		eigenvalue = power_iter(A);
		std::cout << " s3rmt3m3.mtx: eigenvalue from direct power iteration: " << '\t' << eigenvalue << std::endl;

	}
	else if (implem_choice == 2) { // This means applying the power iteration to the Lanczos matrix obtained from
		// s3rmt3m3.mtx
		A = A.convert(argv[2]);
		m=atoi(argv[3]); // convert from char to int
		threshold = atof(argv[4]); // convert from char to double
		eigenvalue = lanczos(A, m, threshold);
		std::cout << " s3rmt3m3.mtx: eigenvalue from Lanczos: " << '\t' << eigenvalue << std::endl;
	}


	return 0;



}
