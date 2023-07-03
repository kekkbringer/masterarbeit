#pragma once

#include <Eigen/Core>

//Eigen::VectorXcd cphf(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, const Eigen::VectorXcd& F);
std::vector<Eigen::VectorXcd> cphf(Eigen::MatrixXcd A, Eigen::MatrixXcd B, std::vector<Eigen::VectorXcd> F);
void tdhf(Eigen::MatrixXcd A, Eigen::MatrixXcd B);
