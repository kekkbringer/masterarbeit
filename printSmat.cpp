#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <vector>
#include "read_herm.hpp"
#include "read_spinor.hpp"

int main() {
	//const auto smat = readHerm("hmat");
	std::vector<double> epsilon;
	const auto spinor = readSpinor(epsilon);
	//std::cout << std::fixed << std::setprecision(8) << spinor << "\n\n";
	std::cout << "real:\n";
	for (int i=0; i<spinor.rows(); i++) {
		for (int j=0; j<spinor.cols(); j++) {
			std::cout << std::fixed << std::setprecision(8) << spinor(j, i).real() << "\n";
		}
	}
	std::cout << "\n\n\nimag:\n";
	for (int i=0; i<spinor.rows(); i++) {
		for (int j=0; j<spinor.cols(); j++) {
			std::cout << std::fixed << std::setprecision(8) << spinor(j, i).imag() << "\n";
		}
	}
}
