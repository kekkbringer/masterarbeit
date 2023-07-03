#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <math.h>

Eigen::MatrixXcd exp(Eigen::MatrixXcd mat) {
	Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(mat.rows(), mat.cols());
	Eigen::MatrixXcd tmp = Eigen::MatrixXcd::Identity(mat.rows(), mat.cols());

	for (int k=0; k<25; k++) {
		res += 1.0 / tgamma(k+1) * tmp;
		tmp *= mat;
	}

	return res;
}

int main() {
	srand((unsigned int) time(0));
	// generate random Matrix
	Eigen::MatrixXcd m = Eigen::MatrixXcd::Random(3, 3);
	// make it hermitian
	auto kek = m + m.adjoint();
	std::cout << "A random, hermitian Matrix:\n" << kek << "\n\n";
	const auto u = exp(1.0j * kek);
	std::cout << "\nexp(i*mat):\n" << u << "\n\n";
	std::cout << "\nU.adjoint * U:\n" << u.adjoint() * u << "\n\n";

}
