#pragma once

#include <Eigen/Core>

Eigen::VectorXcd& eritrans(const Eigen::MatrixXcd &spinorCAO, const int nocc, const int nvirt,
		Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, const Eigen::MatrixXcd &spinor);
