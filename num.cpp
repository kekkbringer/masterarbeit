#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "read_spinor.hpp"

int main() {
	std::vector<double> eref;
	const auto refSpinor = readSpinorDebug(eref, "spinor");
	std::vector<double> eplus;
	const auto plusSpinor = readSpinorDebug(eplus, "../x0plus/spinor");
	std::vector<double> eminus;
	const auto minusSpinor = readSpinorDebug(eminus, "../x0minus/spinor");

	//std::cout << "ref:\n" << refSpinor << "\n";
	//std::cout << "plus:\n" << plusSpinor << "\n";

	// diff plus - ref
	const auto pk = (plusSpinor + refSpinor) / 1e-4;
	std::cout << "plus - ref:\n" << pk << "\n";

	// diff ref - minus
	const auto mk = (refSpinor - minusSpinor) / 1e-4;
	std::cout << "ref - minus:\n" << mk << "\n";

	// diff plus - minus
	const auto pm = (plusSpinor + minusSpinor) / 2e-4;
	std::cout << "plus - minus:\n" << pm << "\n";

	
}
