#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <Eigen/Core>

#include "read_spinor.hpp"

#define DEBUG_MOVEGENERATOR 0
#define IFDBG if constexpr (DEBUG_MOVEGENERATOR)

Eigen::MatrixXcd readSpinor(std::vector<double> &epsilon) {
	// return variable
	Eigen::MatrixXcd s;

	// string for file output, line by line, no eol character
	std::string line;

	// read number of AOs from control file
	// e.g.
	// $rundimensions                                                                                                          
	//    natoms=3                                                                                                             
	//    nbf(CAO)=19                                                                                                          
	//    nbf(AO)=18
	//
	//std::cout << "reading in number of AOs from control file under keyword '$rundimensions'\n";
	std::ifstream control("control");
	if (control.fail()) { std::cout << "\nWARNING: could not read control file\n\n"; return s; }
	int nAO = 0;
	while (getline(control, line)) {
		if (line.substr(0, 14) == "$rundimensions") {
			while (getline(control, line)) {
				std::istringstream iss(line);
				std::string kek;
				while (iss >> kek) {
					if (kek.substr(0, 8) == "nbf(AO)=") {
						nAO = stoi(kek.substr(8));
						break;
					}
				}
			}
		}
	}
	//std::cout << " number of AOs: " << nAO << "\n";
	control.close();


	// return variable
	s.resize(2*nAO, 2*nAO);

	//std::cout << "\n";

	// reading from file spinor.r
	//std::cout << "reading spinor.r ...\n";
	std::ifstream spinorRe("spinor.r");
	if (spinorRe.fail()) {
		std::cout << "\nWARNING: could not read file spinor.r ...\n"
				<< " does it exist?\n do you have permission to read?\n";
		return s;
	}

	// read first line of spinor.r to get format
	// e.g.: $spinor_real    scfconv=7     format(4d20.14)
	getline(spinorRe, line);
	IFDBG std::cout << " first line in spinor.r:\n      " << line << "\n";
	
	// split line by whitespace/tabs/whatever
	std::istringstream iss(line);
	std::string fmt;

	// expecting '$spinor_real'
	iss >> fmt;
	if (fmt != "$spinor_real") {
		std::cout << "\nWARNING: expected flag '$spinor_real' but got: " << fmt << "\n"
				<< " maybe this is the wrong file?\n continuing anyway...\n\n";
	}

	// expecting scfconv=..., doesnt matter for me
	iss >> fmt;

	// expecting format(...), important, will remain in variable 'fmt'
	iss >> fmt;
	IFDBG std::cout << " format: " << fmt << "\n";

	int coeffLength  = 0; // length of coefficients in spinor file
	int coeffPerLine = 0; // number of coefficients given per line
	
	if (fmt.substr(0,6) == "format") {
		// stuff inside parentheses
		std::string format = fmt.substr(7, fmt.length()-2-6);
		IFDBG std::cout << " actual format string: " << format << "\n";
		// iterate over format string
		int i = 0;
		// get coeffPerLine
		while (format[i] != 'f' && format[i] != 'd' && format[i] != 'F' && format[i] != 'D') {
			coeffPerLine *= 10;
			coeffPerLine += format[i] - '0';
			i++;
		}
		IFDBG std::cout << " coefficients per line: " << coeffPerLine << "\n";
		i++; // f or d or F or D, doesnt matter here, we will always use doubles
		// get coeffLength
		while (i < format.length() && format[i] != '.') {
			coeffLength *= 10;
			coeffLength += format[i] - '0';
			i++;
		}
		IFDBG std::cout << " coefficient length: " << coeffLength << "\n";
	} else {
		std::cout << "\nWARNING: could not read format of coefficients: " << fmt << "\n\n";
		return s;
	}

	// reading actualy coefficients
	int i = 0;
	while (getline(spinorRe, line)) {
		if (line[0] == '#') continue;
		if (line.substr(0, 4) == "$end") break;
		// line currently holds spinor irrep, number, eigenvalue, etc...
		std::istringstream energy(line);
		std::string word;
		energy >> word; // number
		energy >> word; // irrep
		energy >> word; // 'eigenvalue=...'
		//std::cout << "orbital energy: " << stod(word.substr(11)) << "\n";
		epsilon.push_back(stod(word.substr(11)));
		for (int mu=0; mu<=2*nAO/coeffPerLine; mu++) {
			getline(spinorRe, line);
			std::replace(line.begin(), line.end(), 'D', 'e'); // convert fortran style to c++
			std::replace(line.begin(), line.end(), 'd', 'e'); // convert fortran style to c++
			//std::cout << mu << ": " << line << "\n";
			//s(i, coeffPerLine*mu + 0) = std::complex<double>(stod(line.substr(0, coeffLength)), 0); }
			//s(i, coeffPerLine*mu + 1) = std::complex<double>(stod(line.substr(coeffLength * 1, coeffLength)), 0); }
			//s(i, coeffPerLine*mu + 2) = std::complex<double>(stod(line.substr(coeffLength * 2, coeffLength)), 0); }
			//s(i, coeffPerLine*mu + 3) = std::complex<double>(stod(line.substr(coeffLength * 3, coeffLength)), 0); }
			for (int k=0; k<coeffPerLine && 2*nAO>coeffPerLine*mu+k; k++) {
				//std::cout << "mu = " << mu << "   k = " << k << "     index = " << (coeffPerLine*mu+k) << "\n";
				//std::cout << "		" << line.substr(coeffLength*k, coeffLength) << "\n";
				s(i, coeffPerLine*mu + k) = std::complex<double>(stod(line.substr(coeffLength*k, coeffLength)), 0);
			}
			if (2*nAO == coeffPerLine * (mu+1)) break;
		}
		i++;
	}
	spinorRe.close();
	//std::cout << "spinor.r read succesfully...\n";






	// reading from file spinor.i
	//std::cout << "reading spinor.i ...\n";
	std::ifstream spinorIm("spinor.i");
	if (spinorIm.fail()) {
		std::cout << "\nWARNING: could not read file spinor.i ...\n"
				<< " does it exist?\n do you have permission to read?\n";
		return s;
	}

	// read first line of spinor.i to get format
	// e.g.: $spinor_imag    scfconv=7     format(4d20.14)
	getline(spinorIm, line);
	IFDBG std::cout << " first line in spinor.i:\n      " << line << "\n";
	
	// split line by whitespace/tabs/whatever
	std::istringstream issi(line);

	// expecting '$spinor_imag'
	issi >> fmt;
	if (fmt != "$spinor_imag") {
		std::cout << "\nWARNING: expected flag '$spinor_imag' but got: " << fmt << "\n"
				<< " maybe this is the wrong file?\n continuing anyway...\n\n";
	}

	// expecting scfconv=..., doesnt matter for me
	issi >> fmt;

	// expecting format(...), important, will remain in variable 'fmt'
	issi >> fmt;
	IFDBG std::cout << " format: " << fmt << "\n";

	coeffLength  = 0; // length of coefficients in spinor file
	coeffPerLine = 0; // number of coefficients given per line
	
	if (fmt.substr(0,6) == "format") {
		// stuff inside parentheses
		std::string format = fmt.substr(7, fmt.length()-2-6);
		IFDBG std::cout << " actual format string: " << format << "\n";
		// iterate over format string
		int i = 0;
		// get coeffPerLine
		while (format[i] != 'f' && format[i] != 'd' && format[i] != 'F' && format[i] != 'D') {
			coeffPerLine *= 10;
			coeffPerLine += format[i] - '0';
			i++;
		}
		IFDBG std::cout << " coefficients per line: " << coeffPerLine << "\n";
		i++; // f or d or F or D, doesnt matter here, we will always use doubles
		// get coeffLength
		while (i < format.length() && format[i] != '.') {
			coeffLength *= 10;
			coeffLength += format[i] - '0';
			i++;
		}
		IFDBG std::cout << " coefficient length: " << coeffLength << "\n";
	} else {
		std::cout << "\nWARNING: could not read format of coefficients: " << fmt << "\n\n";
		return s;
	}

	// reading actualy coefficients
	i = 0;
	while (getline(spinorIm, line)) {
		if (line[0] == '#') continue;
		if (line.substr(0, 4) == "$end") break;
		// line currently holds spinor irrep, number, eigenvalue, etc...
		//std::cout << "\nspinor info: " << line << "\n";
		for (int mu=0; mu<=2*nAO/coeffPerLine; mu++) {
			getline(spinorIm, line);
			std::replace(line.begin(), line.end(), 'D', 'e'); // convert fortran style to c++
			std::replace(line.begin(), line.end(), 'd', 'e'); // convert fortran style to c++
			for (int k=0; k<coeffPerLine && 2*nAO>coeffPerLine*mu+k; k++) {
				s(i, coeffPerLine*mu + k) += std::complex<double>(0, stod(line.substr(coeffLength*k, coeffLength)));
			}
			if (2*nAO == coeffPerLine * (mu+1)) break;
		}
		i++;
	}
	spinorIm.close();
	IFDBG std::cout << "spinor.i read succesfully...\n\n\n";

	return s.transpose();
}
