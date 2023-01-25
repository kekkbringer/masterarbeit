#pragma once

#include <Eigen/Core>
#include "read_fourcenter.hpp"

Eigen::VectorXcd& calcStabmat(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B);
