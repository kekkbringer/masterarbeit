#pragma once

#include <string>
#include <Eigen/Core>

Eigen::MatrixXcd readHerm(std::string f);
//Eigen::MatrixXcd readSymm(std::string f);
Eigen::MatrixXcd readHerm(std::string re, std::string im);
Eigen::MatrixXcd readSymm(std::string re, std::string im);
void splitDipole();
void deleteDipole();
