#include <string>
#include <iostream>
#include <complex>
#include <Eigen/Core>
#include <fstream>
#include <iomanip>
#include <cstdio>

#include "misc.hpp"
#include "read_spinor.hpp"
#include "read_herm.hpp"

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
		if (a=="na") return 11;
		if (a=="mg") return 12;
		if (a=="al") return 13;
		if (a=="si") return 14;
		if (a=="p" ) return 15;
		if (a=="s" ) return 16;
		if (a=="cl") return 17;
		if (a=="ar") return 18;
		return 0;
}

double massOf(std::string a) {
		if (a=="h" ) return 1.008;
		if (a=="he") return 4.0026;
		if (a=="li") return 6.94;
		if (a=="be") return 9.0122;
		if (a=="b" ) return 10.81;
		if (a=="c" ) return 12.011;
		if (a=="n" ) return 14.007;
		if (a=="o" ) return 15.999;
		if (a=="f" ) return 18.998;
		if (a=="ne") return 20.180;
		if (a=="na") return 22.990;
		if (a=="mg") return 24.305;
		if (a=="al") return 26.982;
		if (a=="si") return 28.085;
		if (a=="p" ) return 30.974;
		if (a=="s" ) return 32.06;
		if (a=="cl") return 35.45;
		if (a=="ar") return 39.948;
		return 0.0;
}

void info(int& atomNumber, int& noccupied, int& nvirtual, double& bfieldx, double& bfieldy, double& bfieldz, double& bfieldnorm) {
	std::vector<double> epsilon;
	const auto spinor = readSpinor(epsilon);
	const auto spinorSize = spinor.rows();
	/*****************************************************************************
	 *                          get info about molecule                          *
	 *****************************************************************************
	 *  - read coord file to get order of atoms
	 *  - read basis file to get number of basis functions that are located on
	 *  		each atom
	 */

	// reading coord file to know which atoms occur in which order
	std::ifstream coord("coord");
	if (coord.fail()) std::cout << "\nWARNING: could not read coord file!\n";
	std::string line;
	getline(coord, line); // $coord
	std::vector<std::string> atoms;
	std::vector<double> coordx;
	std::vector<double> coordy;
	std::vector<double> coordz;
	while(getline(coord, line)) {
		if (line[0] == '#') continue;
		if (line[0] == '$') break;
		std::string word;
		std::istringstream iss(line);
		//while(iss >> word) {} // now word contains the elemental symbol
		iss >> word;
		coordx.push_back(std::stod(word));
		iss >> word;
		coordy.push_back(std::stod(word));
		iss >> word;
		coordz.push_back(std::stod(word));
		//std::cout << "elemental symbol: " << word << "\n";
		iss >> word;
		atoms.push_back(word);
	}
	coord.close();
	const int atomNum = atoms.size();


	// get number of basis functions of each atom
	std::ifstream basis("basis");
	if (basis.fail()) std::cout << "\nWARNING: could not read basis file!\n";
	//std::cout << "      minor warning: I'm reading from comments in basis file, e.g. '# li  (7s3p) / [3s2p]    {511/21}'\n";
	//std::cout << "      minor warning: only s-, p- and d-type orbitals are supported yet!\n";

	std::vector<int> s(atomNum);
	std::vector<int> p(atomNum);
	std::vector<int> d(atomNum);
	std::vector<int> f(atomNum);
	int sSize = 0;
	int pSize = 0;
	int dSize = 0;
	int fSize = 0;
	//std::cout << " looking for comments\n";
	while (getline(basis, line)) {
		if (line[0] == '#') {
			//std::cout << line << "\n";
			std::istringstream bss(line);
			std::string word, a;
			bss >> word; // '#'
			bss >> a; // 'li'
			//std::cout << "reading basis info for: " << a << "\n";
			while (bss >> word) {
				if (word[0] == '[') {
					//std::cout << " yesh: " << word << "\n";
					// check for atom type
					for (int atom=0; atom<atomNum; atom++) {
						if (atoms[atom] == a) { // hit
							// iterate over expression like [3s2p]
							int i = 1;
							s.push_back(0);
							for (; i<word.length()-1; i++) { // s
								if (word[i] == 's') break;
								s[atom] *= 10;
								s[atom] += word[i] - '0';
							}
							sSize += s[atom];
							i++;
							p.push_back(0);
							for (; i<word.length()-1; i++) { // p
								if (word[i] == 'p') break;
								p[atom] *= 10;
								p[atom] += word[i] - '0';
							}
							pSize += p[atom];
							i++;
							d.push_back(0);
							for (; i<word.length()-1; i++) { // d
								if (word[i] == 'd') break;
								d[atom] *= 10;
								d[atom] += word[i] - '0';
							}
							dSize += d[atom];
							i++;
							f.push_back(0);
							for (; i<word.length()-1; i++) { // d
								if (word[i] == 'f') break;
								f[atom] *= 10;
								f[atom] += word[i] - '0';
							}
							fSize += f[atom];
						}
					}
				}
			}
		}
	}
	std::vector<int> basisFunctionLength(atomNum);
	int matrixSize = 0;
	for (int i=0; i<atomNum; i++) {
		basisFunctionLength[i] = s[i] + 3*p[i] + 5*d[i] + 7*f[i];
		matrixSize += s[i] + 3*p[i] + 5*d[i] + 7*f[i];
	}
	//std::cout << "matrix size: " << matrixSize << "\n\n";



	//std::cout << "      number of atoms: " << atomNum << "\n";
	//int lower = 0;
	//int upper = 0;
	//for (int i=0; i<atomNum; i++)  {
	//	//std::cout << atoms[i] << "   " << basisFunctionLength[i] << "\n";
	//	if (i<nuc) lower += basisFunctionLength[i];
	//	if (i<=nuc) upper += basisFunctionLength[i];
	//}
	//std::cout << "lower = " << lower << "\n";
	//std::cout << "upper = " << upper << "\n";


	/**********************************************************************
	 *                         read control file                          *
	 **********************************************************************
	 * 
	 * current infos taken from control file:
	 * 	- number of occupied spinors
	 */

	// looking for occupied spinors in control file
	//std::cout << "\n:: reading control file for some additional infos...\n" << std::flush;
	int nocc = 0;
	double Bx = 0.0;
	double By = 0.0;
	double Bz = 0.0;
	double Bnorm = 0.0;
	std::ifstream control("control");
	if (control.fail()) {std::cout << "\nWARNING: cant read control file!\n";}
	while (getline(control, line)) {
		if (line == "$spinor shells") {
			getline(control, line);
			//std::cout << line << "\n"; // e.g. a    1-4                     ( 1 )
			std::istringstream iss(line);
			std::string word;
			iss >> word; // 'a'
			iss >> word;
			//std::cout << "\noccupied spinors: " << word << "\n";
			//std::cout << "  '-' at position " << word.find("-") << "\n";
			//std::cout << "  first number:  " << word.substr(0, word.find("-")) << "\n";
			//std::cout << "  second number: " << word.substr(word.find("-")+1) << "\n";
			nocc = stoi(word.substr(word.find("-")+1)) - stoi(word.substr(0, word.find("-"))) + 1;
			//IFDBG std::cout << "  first " << nocc << " spinors are occupied\n";
		} else if (line == "$magnetic field") {
			getline(control, line);
			std::istringstream iss(line);
			std::string word;
			iss >> word;
			Bx = std::stod(word);
			iss >> word;
			By = std::stod(word);
			iss >> word;
			Bz = std::stod(word);
			iss >> word;
			Bnorm = std::stod(word);
		}
	}
	//std::cout << "   first " << nocc << " spinors are occupied\n";
	//std::cout << "   magnetic field strength: " << Bnorm << "\n";
	//std::cout << "      direction:             " << Bx << "  " << By << "  " << Bz << "\n";
	// normalize B
	const double lambda = Bnorm / sqrt(Bx*Bx + By*By + Bz*Bz);
	Bx *= lambda;
	By *= lambda;
	Bz *= lambda;
	//std::cout << "scaled B vector: " << Bx << "  " << By << "  " << Bz << "\n";

	// end of control file section
	control.close();


	//void info(int& atomNumber, int& nocc, int& nvirt, double& bfieldx, double& bfieldy, double& bfieldz, double& bfieldnorm) {
	atomNumber = atomNum;
	noccupied = nocc;
	nvirtual = spinorSize - nocc;
	bfieldx = Bx;
	bfieldy = By;
	bfieldz = Bz;
	bfieldnorm = Bnorm;
}


Eigen::VectorXcd readVector(std::string file) {
	//std::cout << "READVECTOR CALLED WITH ARG " << file << std::flush;
	std::ifstream re(file + ".r");
	std::ifstream im(file + ".i");

	std::string lineRe;
	std::string lineIm;
	getline(re, lineRe);
	getline(im, lineIm);
	const int size = std::stoi(lineRe);

	Eigen::VectorXcd vec(size);

	for (int i=0; i<size; i++) {
		getline(re, lineRe);
		getline(im, lineIm);
		vec(i) = std::complex<double>(std::stod(lineRe), std::stod(lineIm));
	}

	re.close();
	im.close();

	//std::cout << "READVECTOR DONE.\n" << std::flush;
	return vec;
}


void saveVector(Eigen::VectorXcd vec, std::string name) {
	std::ofstream re(name + ".r");
	std::ofstream im(name + ".i");

	re << vec.rows() << "\n";
	im << vec.rows() << "\n";
	for (int i=0; i<vec.rows(); i++) {
		re << std::setprecision(12) << vec(i).real() << "\n";
		im << std::setprecision(12) << vec(i).imag() << "\n";
	}

	re.close();
	im.close();
}


void splitBraKet(int atomNum) {
	std::ifstream re("sbraket.r");
	std::ifstream im("sbraket.i");
	if (re.fail()) std::cout << "coulnd not read braket file!\n";
	if (im.fail()) std::cout << "coulnd not read braket file!\n";
	std::string lineRe;
	std::string lineIm;
	std::string word;

	getline(re, lineRe);
	getline(im, lineIm);
	std::istringstream iss(lineRe);
	iss >> word;
	const int size = std::stoi(word);
	
	const std::string cartDict[] = {"x", "y", "z"};

	for (int I=0; I<atomNum; I++) {
		for (int J=0; J<atomNum; J++) {
			for (int alpha=0; alpha<3; alpha++) {
				for (int beta=0; beta<3; beta++) {
					//std::cout << I << cartDict[alpha] << " | " << J << cartDict[beta] << "\n";
					std::ofstream outRe("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta] + ".r");
					std::ofstream outIm("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta] + ".i");
					outRe << size << "\n";
					outIm << size << "\n";
					for (int s=0; s<size*size; s++) {
						getline(re, lineRe);
						getline(im, lineIm);
						outRe << std::setprecision(12) << lineRe << "\n";
						outIm << std::setprecision(12) << lineIm << "\n";
					}
				}
			}
		}
	}

	re.close();
	im.close();
}

void splitBraKetold(int atomNum) {
	std::ifstream re("sbraket.r");
	std::ifstream im("sbraket.i");
	if (re.fail()) std::cout << "coulnd not read braket file!\n";
	if (im.fail()) std::cout << "coulnd not read braket file!\n";
	std::string lineRe;
	std::string lineIm;

	getline(re, lineRe);
	getline(im, lineIm);
	const int size = std::stoi(lineRe);

	const std::string cartDict[] = {"x", "y", "z"};

	for (int I=0; I<atomNum; I++) {
		for (int J=0; J<atomNum; J++) {
			for (int alpha=0; alpha<3; alpha++) {
				for (int beta=0; beta<3; beta++) {
					//std::cout << I << cartDict[alpha] << " | " << J << cartDict[beta] << "\n";
					std::ofstream outRe("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta] + ".r");
					std::ofstream outIm("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta] + ".i");
					outRe << size << "\n";
					outIm << size << "\n";
					for (int s=0; s<size; s++) {
						getline(re, lineRe);
						getline(im, lineIm);
						outRe << std::setprecision(12) << lineRe << "\n";
						outIm << std::setprecision(12) << lineIm << "\n";
					}
				}
			}
		}
	}

	re.close();
	im.close();
}

Eigen::MatrixXcd readNumSpinor(std::string filename) {
	std::ifstream file(filename);
	if (file.fail()) std::cout << "COULD NOT READ FILE!!!\n\n";
	std::string line;
	std::string word;
	double re;
	double im;
	getline(file, line);
	std::istringstream iss(line);
	iss >> word;
	const int spinorSize = std::stoi(word);
	//std::cout << "spinorSize = " << spinorSize << "\n";
	Eigen::MatrixXcd res(spinorSize, spinorSize);
	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			getline(file, line);
			//std::cout << "line: " << line << "\n";
			unsigned first = line.find("(");
			unsigned middle = line.find(",");
			unsigned last = line.find(")");

			word = line.substr(first+1, middle-first);
			//std::cout << "   " << word << "\n";
			re = std::stod(word);

			word = line.substr(middle+1, last-middle);
			//std::cout << "   " << word << "\n";
			im = std::stod(word);

			res(j, i) = std::complex<double>(re, im);
		}
	}
	return res;
}

Eigen::MatrixXcd readMatrix(std::string filename) {
	//std::cout << "READMATRIX CALLED WITH ARG " << filename << std::flush;
	using namespace std::complex_literals;

	std::ifstream refile(filename + ".r");
	std::ifstream imfile(filename + ".i");
	if (refile.fail() or imfile.fail()) std::cout << "\n\n\nCOUND NOT READ FILE!!!\n\n\n";

	std::string line, word, reline, imline;
	double re, im;

	getline(refile, line);
	getline(imfile, line);
	std::istringstream iss(line);
	iss >> word; // nlambda
	const int nlambda = std::stoi(word);
	iss >> word; // natom * 3 oder so
	
	Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(nlambda, nlambda);

	for (int i=0; i<nlambda; i++) {
		for (int j=0; j<nlambda; j++) {
			getline(refile, reline);
			getline(imfile, imline);
			//res(i, j) = std::stod(reline) + std::stod(imline) * 1.0i;
			try {res(i, j) += std::stod(reline);} catch (...) {}
			try {res(i, j) += 1.0i * std::stod(imline);} catch (...) {}
		}
	}

	//std::cout << "READMATRIX DONE.\n\n" << std::flush;
	return res.conjugate();
}

Eigen::MatrixXcd readMatrixTransform(std::string filename) {
	// read berrytrans.r for transformation matrix
	std::ifstream transfile("berrytrans.r");
	if (transfile.fail()) throw std::runtime_error("could not find transformation matrix!\n");
	std::string line;
	getline(transfile, line);
	const int dima = std::stoi(line);
	getline(transfile, line);
	const int dimb = std::stoi(line);
	Eigen::MatrixXcd transMat = Eigen::MatrixXcd::Zero(dima, dimb);
	for (int i=0; i<dima; i++) {
		for (int j=0; j<dimb; j++) {
			getline(transfile, line);
			transMat(i, j) += std::stod(line);
		}
	}
	transfile.close();
	//std::cout << transMat.real() << "\n";
	

	// read berryswap.r for reordering matrix
	std::ifstream swapfile("berryswap.r");
	if (swapfile.fail()) throw std::runtime_error("could not find swap matrix!\n");
	getline(swapfile, line);
	const int dims = std::stoi(line);
	Eigen::MatrixXcd swapMat = Eigen::MatrixXcd::Zero(dims, dims);
	for (int i=0; i<dims; i++) {
		for (int j=0; j<dims; j++) {
			getline(swapfile, line);
			swapMat(i, j) += std::stod(line);
		}
	}
	swapfile.close();
	//std::cout << swapMat.real() << "\n";

	//std::cout << "READMATRIX CALLED WITH ARG " << filename << std::flush;
	using namespace std::complex_literals;

	std::ifstream refile(filename + ".r");
	std::ifstream imfile(filename + ".i");
	if (refile.fail() or imfile.fail()) std::cout << "\n\n\nCOUND NOT READ FILE!!!\n\n\n";

	std::string word, reline, imline;
	double re, im;

	getline(refile, line);
	getline(imfile, line);
	std::istringstream iss(line);
	iss >> word; // nlambda
	const int nlambda = std::stoi(word);
	iss >> word; // natom * 3 oder so
	
	Eigen::MatrixXcd res = Eigen::MatrixXcd::Zero(nlambda, nlambda);

	for (int i=0; i<nlambda; i++) {
		for (int j=0; j<nlambda; j++) {
			getline(refile, reline);
			getline(imfile, imline);
			//res(i, j) = std::stod(reline) + std::stod(imline) * 1.0i;
			try {res(i, j) += std::stod(reline);} catch (...) {}
			try {res(i, j) += 1.0i * std::stod(imline);} catch (...) {}
		}
	}

	//std::cout << "READMATRIX DONE.\n\n" << std::flush;
	return transMat.transpose() * swapMat.transpose() * res.conjugate() * swapMat * transMat;
}

Eigen::MatrixXcd readFDEBUG(std::string filename) {
	// read berrytrans.r for transformation matrix
	std::ifstream transfile("berrytrans.r");
	if (transfile.fail()) throw std::runtime_error("could not find transformation matrix!\n");
	std::string line;
	getline(transfile, line);
	const int dima = std::stoi(line);
	getline(transfile, line);
	const int dimb = std::stoi(line);
	Eigen::MatrixXcd transMat = Eigen::MatrixXcd::Zero(dima, dimb);
	for (int i=0; i<dima; i++) {
		for (int j=0; j<dimb; j++) {
			getline(transfile, line);
			transMat(i, j) += std::stod(line);
		}
	}
	transfile.close();
	//std::cout << transMat.real() << "\n";
	

	// read berryswap.r for reordering matrix
	std::ifstream swapfile("berryswap.r");
	if (swapfile.fail()) throw std::runtime_error("could not find swap matrix!\n");
	getline(swapfile, line);
	const int dims = std::stoi(line);
	Eigen::MatrixXcd swapMat = Eigen::MatrixXcd::Zero(dims, dims);
	for (int i=0; i<dims; i++) {
		for (int j=0; j<dims; j++) {
			getline(swapfile, line);
			swapMat(i, j) += std::stod(line);
		}
	}
	swapfile.close();
	//std::cout << swapMat.real() << "\n";

	//std::cout << "READMATRIX CALLED WITH ARG " << filename << std::flush;
	using namespace std::complex_literals;

	auto res = readHerm(filename + ".r", filename + ".i");

	std::cout << "mat dim:   " << res.cols() << " x " << res.rows() << "\n";
	std::cout << "swap dim:  " << swapMat.cols() << " x " << swapMat.rows() << "\n";
	std::cout << "trans dim: " << transMat.cols() << " x " << transMat.rows() << "\n";
	std::cout << std::flush;
	
	//std::cout << "READMATRIX DONE.\n\n" << std::flush;
	return transMat.transpose() * swapMat.transpose() * res.conjugate() * swapMat * transMat;
}

void deleteTmpFiles(const int natom) {
	std::remove("berryswap.r");
	std::remove("berryswap.i");
	std::remove("berrytrans.r");
	std::remove("berrytrans.i");
	for (int I=0; I<natom; I++) {
		std::remove(("b0ai" + std::to_string(I) + "_0.r").c_str());
		std::remove(("b0ai" + std::to_string(I) + "_0.i").c_str());
		std::remove(("b0ai" + std::to_string(I) + "_1.r").c_str());
		std::remove(("b0ai" + std::to_string(I) + "_1.i").c_str());
		std::remove(("b0ai" + std::to_string(I) + "_2.r").c_str());
		std::remove(("b0ai" + std::to_string(I) + "_2.i").c_str());

		std::remove(("b" + std::to_string(I) + "x.r").c_str());
		std::remove(("b" + std::to_string(I) + "x.i").c_str());
		std::remove(("b" + std::to_string(I) + "y.r").c_str());
		std::remove(("b" + std::to_string(I) + "y.i").c_str());
		std::remove(("b" + std::to_string(I) + "z.r").c_str());
		std::remove(("b" + std::to_string(I) + "z.i").c_str());

		std::remove(("h" + std::to_string(I) + "x.r").c_str());
		std::remove(("h" + std::to_string(I) + "x.i").c_str());
		std::remove(("h" + std::to_string(I) + "y.r").c_str());
		std::remove(("h" + std::to_string(I) + "y.i").c_str());
		std::remove(("h" + std::to_string(I) + "z.r").c_str());
		std::remove(("h" + std::to_string(I) + "z.i").c_str());

		std::remove(("k" + std::to_string(I) + "x.r").c_str());
		std::remove(("k" + std::to_string(I) + "x.i").c_str());
		std::remove(("k" + std::to_string(I) + "y.r").c_str());
		std::remove(("k" + std::to_string(I) + "y.i").c_str());
		std::remove(("k" + std::to_string(I) + "z.r").c_str());
		std::remove(("k" + std::to_string(I) + "z.i").c_str());

		std::remove(("snxbraip" + std::to_string(I) + "_0.r").c_str());
		std::remove(("snxbraip" + std::to_string(I) + "_0.i").c_str());
		std::remove(("snxbraip" + std::to_string(I) + "_1.r").c_str());
		std::remove(("snxbraip" + std::to_string(I) + "_1.i").c_str());
		std::remove(("snxbraip" + std::to_string(I) + "_2.r").c_str());
		std::remove(("snxbraip" + std::to_string(I) + "_2.i").c_str());

		std::remove(("snxketip" + std::to_string(I) + "_0.r").c_str());
		std::remove(("snxketip" + std::to_string(I) + "_0.i").c_str());
		std::remove(("snxketip" + std::to_string(I) + "_1.r").c_str());
		std::remove(("snxketip" + std::to_string(I) + "_1.i").c_str());
		std::remove(("snxketip" + std::to_string(I) + "_2.r").c_str());
		std::remove(("snxketip" + std::to_string(I) + "_2.i").c_str());

		for (int J=0; J<natom; J++) {
			std::remove(("bk" + std::to_string(I) + "x" + std::to_string(J) + "x.r").c_str());
			std::remove(("bk" + std::to_string(I) + "y" + std::to_string(J) + "y.r").c_str());
			std::remove(("bk" + std::to_string(I) + "z" + std::to_string(J) + "z.r").c_str());
			std::remove(("bk" + std::to_string(I) + "x" + std::to_string(J) + "y.r").c_str());
			std::remove(("bk" + std::to_string(I) + "y" + std::to_string(J) + "x.r").c_str());
			std::remove(("bk" + std::to_string(I) + "x" + std::to_string(J) + "z.r").c_str());
			std::remove(("bk" + std::to_string(I) + "z" + std::to_string(J) + "x.r").c_str());
			std::remove(("bk" + std::to_string(I) + "y" + std::to_string(J) + "z.r").c_str());
			std::remove(("bk" + std::to_string(I) + "z" + std::to_string(J) + "y.r").c_str());
			std::remove(("bk" + std::to_string(I) + "x" + std::to_string(J) + "x.i").c_str());
			std::remove(("bk" + std::to_string(I) + "y" + std::to_string(J) + "y.i").c_str());
			std::remove(("bk" + std::to_string(I) + "z" + std::to_string(J) + "z.i").c_str());
			std::remove(("bk" + std::to_string(I) + "x" + std::to_string(J) + "y.i").c_str());
			std::remove(("bk" + std::to_string(I) + "y" + std::to_string(J) + "x.i").c_str());
			std::remove(("bk" + std::to_string(I) + "x" + std::to_string(J) + "z.i").c_str());
			std::remove(("bk" + std::to_string(I) + "z" + std::to_string(J) + "x.i").c_str());
			std::remove(("bk" + std::to_string(I) + "y" + std::to_string(J) + "z.i").c_str());
			std::remove(("bk" + std::to_string(I) + "z" + std::to_string(J) + "y.i").c_str());
		}
	}

	std::remove("cfock.r");
	std::remove("cfock.i");
	std::remove("xfock.r");
	std::remove("xfock.i");
	std::remove("edipole.r");
	std::remove("edipole.i");
	std::remove("hmatcao.r");
	std::remove("hmatcao.i");
	//std::remove("smatcao.r");
	//std::remove("smatcao.i");
	std::remove("wcao.r");
	std::remove("wcao.i");
}
