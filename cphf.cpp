#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <typeinfo>
#include <iomanip>

#include "cphf.hpp"

Eigen::VectorXcd cphf(Eigen::MatrixXcd A, Eigen::MatrixXcd B, Eigen::VectorXcd F) {
	Eigen::MatrixXcd kek(2*A.rows(), 2*A.rows());
	//kek << A, B,
	//       B.conjugate(), A.conjugate();
	kek << A.conjugate(), B.conjugate(),
	       B, A;

	Eigen::VectorXcd schmon(2*A.rows());
	schmon << F, F.conjugate();

	Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> dec(kek);
	//Eigen::LLT<Eigen::MatrixXcd> dec(kek);
	auto u = dec.solve(schmon);

	const auto info = dec.info();
	if (info == Eigen::ComputationInfo::Success) {
		std::cout << "      converged!\n";
	} else if (info == Eigen::ComputationInfo::NumericalIssue) {
		std::cout << "      numerical issue! Matrix was not positive definite!\n";
		//std::cout << "        RETURNING 0-VECTOR!\n\n";
		//return Eigen::VectorXcd::Zero(2*A.rows());
	} else if (info == Eigen::ComputationInfo::NoConvergence) {
		std::cout << "      NOT converged!\n";
	} else if (info == Eigen::ComputationInfo::InvalidInput) {
		std::cout << "      unknown error in Cholesky decomposition!\n";
	} 

	const double relativeError = (kek*u - schmon).norm() / u.norm();
	std::cout << "      relative error: " << relativeError << "\n";

	std::cout << "CPHF done\n\n" << std::flush;
	return u;
}

void tdhf(Eigen::MatrixXcd A, Eigen::MatrixXcd B) {
	Eigen::MatrixXcd kek(2*A.rows(), 2*A.rows());
	kek << A, B,
	       B.conjugate(), A.conjugate();

	Eigen::MatrixXcd schmon(2*A.rows(), 2*A.rows());
	schmon << Eigen::MatrixXcd::Identity(A.rows(), A.rows()), Eigen::MatrixXcd::Zero(A.rows(), A.rows()),
		Eigen::MatrixXcd::Zero(A.rows(), A.rows()), -1.0*Eigen::MatrixXcd::Identity(A.rows(), A.rows());

	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ges;
	ges.compute(schmon*kek);

	std::cout << "\n\n:: TDHF excitation energies:\n" << std::setprecision(10);
	std::cout << " relative precision: " << (kek*ges.eigenvectors()).norm() / (ges.eigenvalues().asDiagonal()*schmon*ges.eigenvectors()).norm() << "\n";
	std::cout << ges.eigenvalues() << "\n\n" << std::flush;

}
