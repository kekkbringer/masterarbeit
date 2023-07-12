#include "matrixio.hpp"

#include <fstream>
#include <complex>
#include <iomanip>

#include <iostream>

void saveMatrix(const Eigen::MatrixXcd& mat, const std::string filename) {
	const int rows = mat.rows();
	const int cols = mat.cols();
	
	std::ofstream matfileRe(filename + ".r");
	std::ofstream matfileIm(filename + ".i");

	matfileRe << rows << "\n" << cols << "\n";
	matfileIm << rows << "\n" << cols << "\n";

	for (int i=0; i<rows; i++) {
		for (int j=0; j<cols; j++) {
			matfileRe << std::setprecision(std::numeric_limits<double>::digits10 + 1) << mat.real()(i, j) << "\n";
			matfileIm << std::setprecision(std::numeric_limits<double>::digits10 + 1) << mat.imag()(i, j) << "\n";
		}
	}

	matfileRe.close();
	matfileIm.close();
}

Eigen::MatrixXcd loadMatrix(std::string filename) {
	std::ifstream matfileRe(filename + ".r");
	std::ifstream matfileIm(filename + ".i");

	std::string lineRe, lineIm;

	getline(matfileRe, lineRe);
	getline(matfileIm, lineIm);
	const int rows = stoi(lineRe);
	getline(matfileRe, lineRe);
	getline(matfileIm, lineIm);
	const int cols = stoi(lineRe);

	Eigen::MatrixXcd mat(rows, cols);

	for (int i=0; i<rows; i++) {
		for (int j=0; j<cols; j++) {
			getline(matfileRe, lineRe);
			getline(matfileIm, lineIm);
			mat(i, j) = std::complex<double>(stod(lineRe), stod(lineIm));
		}
	}

	return mat;
}
