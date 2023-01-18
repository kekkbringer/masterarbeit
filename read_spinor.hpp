#pragma once

#include <Eigen/Core>
#include <vector>
#include <string>

Eigen::MatrixXcd readSpinor(std::vector<double> &epsilon);
Eigen::MatrixXcd readMos(std::vector<double> &epsilon);
Eigen::MatrixXcd readSpinorDebug(std::vector<double> &epsilon, std::string file);
