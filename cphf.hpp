#pragma once

#include <Eigen/Core>

Eigen::VectorXcd cphf(Eigen::MatrixXcd A, Eigen::MatrixXcd B, Eigen::VectorXcd F);
void tdhf(Eigen::MatrixXcd A, Eigen::MatrixXcd B);
