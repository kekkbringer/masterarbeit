#pragma once

#include <Eigen/Core>
#include <string>

int chargeOf(std::string a);
double massOf(std::string a);

void info(int& atomNumber, int& noccupied, int& nvirtual, double& bfieldx, double& bfieldy, double& bfieldz, double& bfieldnorm);

Eigen::VectorXcd readVector(std::string file);
Eigen::MatrixXcd readMatrix(std::string file);
Eigen::MatrixXcd readMatrixTransform(std::string file);
Eigen::MatrixXcd readFDEBUG(std::string file);
void saveVector(Eigen::VectorXcd vec, std::string name);
void splitBraKet(int atomNum);
Eigen::MatrixXcd readNumSpinor(std::string filename);
void deleteTmpFiles(const int natom);
Eigen::MatrixXcd readCoulomb(std::string filename);
Eigen::MatrixXd readExchange(std::string filename, int c);
