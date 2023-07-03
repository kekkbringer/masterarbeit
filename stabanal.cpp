#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <complex>

#include "calc_stabmat.hpp"

using namespace std::complex_literals;

int main(int argc, char* argv[]) {
	Eigen::MatrixXcd A;
	Eigen::MatrixXcd B;
	calcStabmat(A, B);

	Eigen::MatrixXcd kek(2*A.rows(), 2*A.rows());
	kek << A, B, B.conjugate(), A.conjugate();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(kek);
	std::cout << "Eigenvalues:\n" << es.eigenvalues() << "\n\n";

	return 0;
}
