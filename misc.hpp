#pragma once

#include <string>
#include <iostream>
#include <complex>
#include <Eigen/Core>

int chargeOf(std::string a) {
		if (a=="h" ) return 1;
		if (a=="he") return 2;
		if (a=="li") return 3;
		if (a=="be") return 4;
		if (a=="b" ) return 5;
		if (a=="c" ) return 6;
		if (a=="n" ) return 7;
		if (a=="o" ) return 8;
		if (a=="f" ) return 9;
		if (a=="ne") return 10;
		return 0;
}
