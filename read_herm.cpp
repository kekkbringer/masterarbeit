#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <Eigen/Core>
#include <cstdio>

void deleteDipole() {
	std::remove("dipolex.r");
	std::remove("dipoley.r");
	std::remove("dipolez.r");
	std::remove("dipolex.i");
	std::remove("dipoley.i");
	std::remove("dipolez.i");
}

Eigen::MatrixXcd readSymm(std::string re, std::string im) {
	Eigen::MatrixXcd mat;

	// reading real part
	//std::cout << "reading file " << f << ".r...\n";
	std::ifstream matRe(re);
	if (matRe.fail()) { std::cout << "WARNING: could not read file " << re << "!\n\n"; return mat; }
	std::ifstream matIm(im);
	if (matIm.fail()) { std::cout << "WARNING: could not read file " << im << "!\n\n"; return mat; }

	std::string lineRe;
	std::string lineIm;
	
	// get "size" of mat (number of elements in upper half)
	getline(matRe, lineRe);
	getline(matIm, lineIm);

	int nAO = -0.5 + sqrt(2.0*std::stoi(lineRe) + 0.25);
	std::cout << "nAO: " << nAO << "\n";

	mat.resize(nAO, nAO);
	for (int i=0; i<nAO; i++) {
		for (int j=0; j<=i; j++) {
			getline(matRe, lineRe);
			getline(matIm, lineIm);
			mat(i, j) = std::complex<double>(stod(lineRe), stod(lineIm));
			mat(j, i) = std::complex<double>(stod(lineRe), stod(lineIm));
		}
	}
	matRe.close();
	matIm.close();

	return mat;
}


//Eigen::MatrixXcd readSymm(std::string f) {
//	Eigen::MatrixXcd mat;
//
//	// string for file output, line by line, no eol character
//	std::string line;
//
//	// read number of AOs from control file
//	// e.g.
//	// $rundimensions                                                                                                          
//	//    natoms=3                                                                                                             
//	//    nbf(CAO)=19                                                                                                          
//	//    nbf(AO)=18
//	//
//	//std::cout << "reading in number of AOs from control file under keyword '$rundimensions'\n";
//	std::ifstream control("control");
//	if (control.fail()) { std::cout << "\nWARNING: could not read control file\n\n"; return mat; }
//	int nAO = 0;
//	while (getline(control, line)) {
//		if (line.substr(0, 14) == "$rundimensions") {
//			while (getline(control, line)) {
//				std::istringstream iss(line);
//				std::string kek;
//				while (iss >> kek) {
//					if (kek.substr(0, 8) == "nbf(AO)=") {
//						nAO = stoi(kek.substr(8));
//						break;
//					}
//				}
//			}
//		}
//	}
//	//std::cout << " number of AOs: " << nAO << "\n";
//	control.close();
//
//	mat.resize(nAO, nAO);
//
//	// reading real part
//	//std::cout << "reading file " << f << ".r...\n";
//	std::ifstream matRe(f + ".r");
//	if (matRe.fail()) { std::cout << "WARNING: could not read file!\n\n"; return mat; }
//	std::ifstream matIm(f + ".i");
//	if (matIm.fail()) { std::cout << "WARNING: could not read file!\n\n"; return mat; }
//
//	std::string lineRe;
//	std::string lineIm;
//	
//	// get "size" of mat (number of elements in upper half)
//	getline(matRe, lineRe);
//	getline(matIm, lineIm);
//
//	int nAO = -0.5 + sqrt(2.0*std::stoi(lineRe) + 0.25);
//	std::cout << "nAO: " << nAO << "\n";
//
//	for (int i=0; i<nAO; i++) {
//		for (int j=0; j<=i; j++) {
//			getline(matRe, lineRe);
//			getline(matIm, lineIm);
//			mat(i, j) = std::complex<double>(stod(lineRe), stod(lineIm));
//			mat(j, i) = std::complex<double>(stod(lineRe), stod(lineIm));
//		}
//	}
//	matRe.close();
//	matIm.close();
//
//	return mat;
//}

Eigen::MatrixXcd readHerm(std::string re, std::string im) {
	Eigen::MatrixXcd mat;

	// reading real part
	//std::cout << "reading file " << f << ".r...\n";
	std::ifstream matRe(re);
	if (matRe.fail()) { std::cout << "WARNING: could not read file " << re << "!\n\n"; return mat; }
	std::ifstream matIm(im);
	if (matIm.fail()) { std::cout << "WARNING: could not read file " << im << "!\n\n"; return mat; }

	std::string lineRe;
	std::string lineIm;
	
	// get "size" of mat (number of elements in upper half)
	getline(matRe, lineRe);
	getline(matIm, lineIm);

	int nAO = -0.5 + sqrt(2.0*std::stoi(lineRe) + 0.25);
	std::cout << "nAO: " << nAO << "\n";

	mat.resize(nAO, nAO);
	for (int i=0; i<nAO; i++) {
		for (int j=0; j<=i; j++) {
			getline(matRe, lineRe);
			getline(matIm, lineIm);
			mat(i, j) = std::complex<double>(stod(lineRe), -stod(lineIm));
			mat(j, i) = std::complex<double>(stod(lineRe),  stod(lineIm));
		}
	}
	matRe.close();
	matIm.close();

	return mat;
}


Eigen::MatrixXcd readHerm(std::string f) {
	Eigen::MatrixXcd mat;

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
	if (control.fail()) { std::cout << "\nWARNING: could not read control file\n\n"; return mat; }
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

	mat.resize(nAO, nAO);

	// reading real part
	//std::cout << "reading file " << f << ".r...\n";
	std::ifstream matRe(f + ".r");
	if (matRe.fail()) { std::cout << "WARNING: could not read file!\n\n"; return mat; }
	std::ifstream matIm(f + ".i");
	if (matIm.fail()) { std::cout << "WARNING: could not read file!\n\n"; return mat; }

	std::string lineRe;
	std::string lineIm;
	
	// get "size" of mat (number of elements in upper half)
	getline(matRe, lineRe);
	getline(matIm, lineIm);
	int comp = nAO*(nAO-1)/2 + nAO;
	//if (stoi(lineRe) != comp) {std::cout << "\nWARNING: dimensions dont match!" << comp << " vs " << lineRe << "\n\n"; return mat;}

	for (int i=0; i<nAO; i++) {
		for (int j=0; j<=i; j++) {
			getline(matRe, lineRe);
			getline(matIm, lineIm);
			mat(i, j) = std::complex<double>(stod(lineRe), -stod(lineIm));
			mat(j, i) = std::complex<double>(stod(lineRe),  stod(lineIm));
		}
	}
	matRe.close();
	matIm.close();

	return mat;
}

void splitDipole() {
	//std::cout << "\nsplitting dipole files...\n";
	std::ifstream dipoleRe("edipole.r");
	if (dipoleRe.fail()) {std::cout << "\nWARNING: could not open edipole.r\n";}
	std::ifstream dipoleIm("edipole.i");
	if (dipoleRe.fail()) {std::cout << "\nWARNING: could not open edipole.i\n";}

	std::string line;

	// real part
	getline(dipoleRe, line);
	int size = stoi(line);
	//std::cout << " size of real part: " << size << "\n";
	// x
	std::ofstream xr("dipolex.r");
	xr << size << "\n";
	for (int i=0; i<size; i++) {
		getline(dipoleRe, line);
		xr << line << "\n";
	}
	// y
	std::ofstream yr("dipoley.r");
	yr << size << "\n";
	for (int i=0; i<size; i++) {
		getline(dipoleRe, line);
		yr << line << "\n";
	}
	// z
	std::ofstream zr("dipolez.r");
	zr << size << "\n";
	for (int i=0; i<size; i++) {
		getline(dipoleRe, line);
		zr << line << "\n";
	}
	dipoleRe.close();


	// imag part
	getline(dipoleIm, line);
	size = stoi(line);
	//std::cout << " size of imag part: " << size << "\n";
	// x
	std::ofstream xi("dipolex.i");
	xi << size << "\n";
	for (int i=0; i<size; i++) {
		getline(dipoleIm, line);
		xi << line << "\n";
	}
	// y
	std::ofstream yi("dipoley.i");
	yi << size << "\n";
	for (int i=0; i<size; i++) {
		getline(dipoleIm, line);
		yi << line << "\n";
	}
	// z
	std::ofstream zi("dipolez.i");
	zi << size << "\n";
	for (int i=0; i<size; i++) {
		getline(dipoleIm, line);
		zi << line << "\n";
	}
	dipoleIm.close();

	xr.close();
	yr.close();
	zr.close();
	xi.close();
	yi.close();
	zi.close();
}
