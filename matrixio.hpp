#pragma once

#include <string>
#include <Eigen/Core>

void saveMatrix(const Eigen::MatrixXcd& mat, const std::string filename);
Eigen::MatrixXcd loadMatrix(std::string filename);
