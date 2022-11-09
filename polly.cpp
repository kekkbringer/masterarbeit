#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <complex>
#include <fstream>
#include <iomanip>

#include "cphf.hpp"
#include "read_spinor.hpp"
#include "read_herm.hpp"
#include "read_fourcenter.hpp"
#include "calc_stabmat.hpp"

#define DEBUG 1
#define IFDBG if constexpr (DEBUG)

#define MAJOR_VERSION 1
#define MINOR_VERSION 0
#define PATCH_VERSION "5-alpha"

int main(int argc, char* argv[]) {	
	/**********************************************************************
	 *                          Checking flags                            *
	 **********************************************************************
	 *
	 * currently suppported flags:
	 * 	-h	help
	 * 	-v	version
	 * 	-tda	switch on Tamm-Dancoff-Approximation
	 */

	bool tda = false;
	bool tdhf_calc = false;
	if (argc > 1) { // flags are present
		if (std::string(argv[1]) == "-h") {
			std::cout << "go ask Dominik.\n" << std::flush;
			return 0;
		} else if (std::string(argv[1]) == "-v") {
			std::cout << "V" << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_VERSION << "\n" << std::flush;
			return 0;
		} else if (std::string(argv[1]) == "-tda") {
			tda = true;
			std::cout << "TDA enabled\n\n";
		} else if (std::string(argv[1]) == "-tdhf") {
			tdhf_calc = true;
			std::cout << "TDHF enabled\n\n";
		}
	}

	//const auto smat = readHerm("smat");
	//std::cout << smat << "\n";
	//return 0;


	/**********************************************************************
	 *                          print banner                              *
	 *********************************************************************/
	std::cout << "      __   __               \n";
	std::cout << "     |__) /  \\ |    |    \\ /\n";
	std::cout << "     |    \\__/ |___ |___  | \n";
	std::cout << "\n" << std::flush;
                        


	/**********************************************************************
	 *                         read control file                          *
	 **********************************************************************
	 * 
	 * current infos taken from control file:
	 * 	- number of occupied spinors
	 */

	bool foundAReal = false;
	bool foundAImag = false;
	bool foundBReal = false;
	bool foundBImag = false;

	std::string aReal = "";
	std::string aImag = "";
	std::string bReal = "";
	std::string bImag = "";

	// looking for occupied spinors in control file
	std::cout << "\n:: reading control file...\n" << std::flush;
	int nocc = 0;
	std::ifstream control("control");
	if (control.fail()) {std::cout << "\nWARNING: cant read control file!\n"; return -1;}
	std::string line;
	while (getline(control, line)) {
		if (line[0] == '#') {
			continue;
		}
		if (line == "$spinor shells") {
			getline(control, line);
			//std::cout << line << "\n"; // e.g. a    1-4                     ( 1 )
			std::istringstream iss(line);
			std::string word;
			iss >> word; // 'a'
			iss >> word;
			nocc = stoi(word.substr(word.find("-")+1)) - stoi(word.substr(0, word.find("-"))) + 1;
		} else if (line.substr(0, 13) == "$polly_a_real") {
			std::istringstream iss(line);
			std::string word;
			iss >> word; // '$polly_a_real'
			iss >> word;
			if (word.substr(0, 5) == "file=") {
				//std::cout << "file sollte " << word.substr(5) << " sein\n";
				aReal = word.substr(5);
				foundAReal = true;
			} else {
				std::cout << "WARNING: found '$polly_a_real' without 'file='!\n";
				std::cout << std::flush;
			}
		} else if (line.substr(0, 13) == "$polly_a_imag") {
			std::istringstream iss(line);
			std::string word;
			iss >> word; // '$polly_a_imag'
			iss >> word;
			if (word.substr(0, 5) == "file=") {
				//std::cout << "file sollte " << word.substr(5) << " sein\n";
				aImag = word.substr(5);
				foundAImag = true;
			} else {
				std::cout << "WARNING: found '$polly_a_imag' without 'file='!\n";
				std::cout << std::flush;
			}
		} else if (line.substr(0, 13) == "$polly_b_real") {
			std::istringstream iss(line);
			std::string word;
			iss >> word; // '$polly_b_real'
			iss >> word;
			if (word.substr(0, 5) == "file=") {
				//std::cout << "file sollte " << word.substr(5) << " sein\n";
				bReal = word.substr(5);
				foundBReal = true;
			} else {
				std::cout << "WARNING: found '$polly_b_real' without 'file='!\n";
				std::cout << std::flush;
			}
		} else if (line.substr(0, 13) == "$polly_b_imag") {
			std::istringstream iss(line);
			std::string word;
			iss >> word; // '$polly_b_imag'
			iss >> word;
			if (word.substr(0, 5) == "file=") {
				//std::cout << "file sollte " << word.substr(5) << " sein\n";
				bImag = word.substr(5);
				foundBImag = true;
			} else {
				std::cout << "WARNING: found '$polly_b_imag' without 'file='!\n";
				std::cout << std::flush;
			}
		}
	}

	//std::cout << "a real: " << aReal << "\n";
	//std::cout << "a imag: " << aImag << "\n";
	//std::cout << "b real: " << bReal << "\n";
	//std::cout << "b imag: " << bImag << "\n";

	// end of control file section
	control.close();

	std::cout << "\n" << std::flush;



	/**********************************************************************
	 *                    Read A and B (if possible)                      *
	 *********************************************************************/
	Eigen::MatrixXcd A;
	Eigen::MatrixXcd B;

	bool calcNeeded = true;
	if (foundAReal and foundAImag and foundBReal and foundBImag) {
		calcNeeded = false;
		std::cout << "Reading Matricies A and B from previous calculation...\n" << std::flush;
	}
	// read A
	if (foundAReal and foundAImag) {
		std::cout << "reading A from files...\n";
		A = readHerm(aReal, aImag);
	}
	// read B
	if (foundBReal and foundBImag) {
		std::cout << "reading B from files...\n";
		B = readSymm(bReal, bImag);
	}



	/**********************************************************************
	 *                          Reading files                             *
	 **********************************************************************
	 *
	 * the following files will be read:
	 * 	- spinor.r
	 * 	- spinor.i
	 * 	- edipole.r
	 * 	- edipole.i
	 */

	std::cout << "\n:: reading files...\n" << std::flush;

	// read dipole files
	std::cout << " reading dipole files...\n" << std::flush;
	splitDipole(); // split file
	const auto dipoleX = readHerm("dipolex");
	const auto dipoleY = readHerm("dipoley");
	const auto dipoleZ = readHerm("dipolez");
	deleteDipole(); // delete split files
	std::cout << "      size of dipoles:   3 x " << dipoleX.rows() << " x " << dipoleX.rows() << "\n";
	std::cout << "      amounts to:        " << 3.0*(sizeof dipoleX)*dipoleX.rows()*dipoleX.rows()/1000000.0 << " MB\n";
	std::cout << "\n" << std::flush;
	
	// read spinor.r and spinor.i files
	std::cout << " reading spinor files...\n" << std::flush;
	std::vector<double> epsilon; // orbital energies spinor
	const auto spinor = readSpinor(epsilon);
	// save spinorsize for matrix size
	const size_t spinorSize = spinor.rows();
	std::cout << "      size of spinor:   " << spinorSize << " x " << spinorSize << "\n";
	std::cout << "      amounts to:       " << (sizeof spinor)*spinorSize*spinorSize/1000000.0 << " MB\n";
	std::cout << "\n" << std::flush;
	const int nvirt = spinorSize-nocc;

	Eigen::VectorXcd Fx(nocc*nvirt);
	Eigen::VectorXcd Fy(nocc*nvirt);
	Eigen::VectorXcd Fz(nocc*nvirt);
		


	/**********************************************************************
	 *                transform matricies into mo-basis                   *
	 **********************************************************************
	 *
	 * 	- dipole moments
	 */

	// transform dipoles from ao basis to spinor basis
	std::cout << "\n:: transforming dipoles to spinor basis...\n" << std::flush;
	Eigen::MatrixXcd dxBig(spinorSize, spinorSize);
	dxBig << dipoleX, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), dipoleX;
	Eigen::MatrixXcd dyBig(spinorSize, spinorSize);
	dyBig << dipoleY, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), dipoleY;
	Eigen::MatrixXcd dzBig(spinorSize, spinorSize);
	dzBig << dipoleZ, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), dipoleZ;

	dxBig = spinor.adjoint() * dxBig * spinor;
	dyBig = spinor.adjoint() * dyBig * spinor;
	dzBig = spinor.adjoint() * dzBig * spinor;

	std::cout << "      size of new dipoles:   3 x " << dxBig.rows() << " x " << dxBig.rows() << "\n";
	std::cout << "      amounts to:            " << 3.0*(sizeof dxBig)*dxBig.rows()*dxBig.rows()/(1000.0*1000.0) << " MB\n";

	// end of transformation section




	/**********************************************************************
	 *                            CPHF section                            *
	 **********************************************************************
	 *
	 * construction of A, B and F matricies for CPHF equaltions:	      */
	//
	//   / A  B \   / U \   / F \
	//   |      | * |   | = |   |
	//   \ B* A*/   \ U*/   \ F*/
	//
	//

	if (calcNeeded) {
		std::cout << "\ncalculating A and B from scratch...\n";
		calcStabmat(A, B);
	} else {
		std::cout << "\nA and B taken from logfiles\n";
	}

	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			const int x = i*nvirt + a - nocc;
			// vector F
			Fx(x) = dxBig(i, a);
			Fy(x) = dyBig(i, a);
			Fz(x) = dzBig(i, a);
		}
	}

	std::cout << "\n" << std::flush;




	/*****************************************************************************
	 *                                  TDHF                                     *
	 *****************************************************************************
	 *
	 * Just for debugging
	 */

	if (tdhf_calc) {
		tdhf(A, B);
		return 0;
	}


	

	/*****************************************************************************
	 *                                  CPHF                                     *
	 *****************************************************************************/
	
	// solving CPHF equations
	std::cout << "\n:: solving CPHF equations...\n" << std::flush;

	std::cout << "      CPHF x...\n";
	const auto ux = cphf(A, B, Fx);
	std::cout << "      CPHF y...\n";
	const auto uy = cphf(A, B, Fy);
	std::cout << "      CPHF z...\n";
	const auto uz = cphf(A, B, Fz);
	
	// end of CPHF section
	
	std::cout << "" << std::flush;


	/**********************************************************************
	 *              calculating and printing final results                *
	 **********************************************************************
	 *
	 * currently only calculating static polarizability
	 */

	Eigen::VectorXcd bigfx(2*nocc*nvirt); bigfx << Fx, Fx.conjugate();
	Eigen::VectorXcd bigfy(2*nocc*nvirt); bigfy << Fy, Fy.conjugate();
	Eigen::VectorXcd bigfz(2*nocc*nvirt); bigfz << Fz, Fz.conjugate();
	
	std::cout << std::fixed;
	std::cout << std::setprecision(8);
	
	std::cout << "\n:: calculating static polarizability...\n" << std::flush;
	const std::complex<double> alphaXX = bigfx.adjoint() * ux;
	const std::complex<double> alphaYY = bigfy.adjoint() * uy;
	const std::complex<double> alphaZZ = bigfz.adjoint() * uz;
	const std::complex<double> alphaXY = bigfy.adjoint() * ux;
	const std::complex<double> alphaXZ = bigfz.adjoint() * ux;
	const std::complex<double> alphaYZ = bigfz.adjoint() * uy;
	const std::complex<double> alphaYX = alphaXY;//bigfx.adjoint() * uy;
	const std::complex<double> alphaZX = alphaXZ;//bigfx.adjoint() * uz;
	const std::complex<double> alphaZY = alphaYZ;//bigfy.adjoint() * uz;

	Eigen::MatrixXd alpha(3, 3);
	alpha << alphaXX.real(), alphaXY.real(), alphaXZ.real(),
	      	 alphaYX.real(), alphaYY.real(), alphaYZ.real(),
		 alphaZX.real(), alphaZY.real(), alphaZZ.real();
	std::cout << "__________________________________________________________________________\n";
	//Eigen::IOFormat def(Eigen::StreamPrecision, 0, "      ", "\n", "             ", "");
	Eigen::IOFormat def(Eigen::StreamPrecision, 0, "      ", "\n", "             ", "");
	std::cout << "      static electronic polarizability:\n" << alpha.format(def) << "\n";
	const auto isotropicAlpha = (alpha(0, 0) + alpha(1, 1) + alpha(2, 2)) / 3.0;
	std::cout << "\n";
	std::cout << "      1/3*Tr(alpha):   " << isotropicAlpha << "\n";
	std::cout << "      Anisotropy:      " << sqrt( 0.5*((alphaXX-alphaYY)*(alphaXX-alphaYY) + (alphaYY-alphaZZ)*(alphaYY-alphaZZ) + (alphaZZ-alphaXX)*(alphaZZ-alphaXX))
						+ 3.0*(alphaXY*alphaXY + alphaYZ*alphaYZ + alphaZX*alphaZX) ).real() << "\n";
	std::cout << "__________________________________________________________________________\n";

	// write result to file 'pola'
	std::ofstream polly("pola");
	polly << std::fixed;
	polly << std::setprecision(8);
	polly << "static electronic polarizability\n";
	polly << alpha.format(def) << "\n";
	polly << "1/3*Tr(alpha):    " << isotropicAlpha << "\n";
	polly << "Anisotropy:       " << sqrt( 0.5*((alphaXX-alphaYY)*(alphaXX-alphaYY) + (alphaYY-alphaZZ)*(alphaYY-alphaZZ) + (alphaZZ-alphaXX)*(alphaZZ-alphaXX))
						+ 3.0*(alphaXY*alphaXY + alphaYZ*alphaYZ + alphaZX*alphaZX) ).real() << "\n";
	



	// return section

	std::cout << "\n\n	polly ended successfully.\n\n" << std::flush;

	return 0;
}
