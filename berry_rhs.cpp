#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <complex>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <stdio.h>

#include "misc.hpp"
#include "cphf.hpp"
#include "read_spinor.hpp"
#include "read_herm.hpp"
#include "read_fourcenter.hpp"
#include "calc_stabmat.hpp"
#include "fci_grad.hpp"

using namespace std::complex_literals;

/*
 * 1 Re(aa)
 * 2 Re(ab+ba)
 * 3 Im(ab+ba)
 * 4 Re(bb)
 * 5 Im(aa)
 * 6 Re(ab-ba)
 * 7 Im(ab-ba)
 * 8 Im(bb)
 *
 * aa = 1 + i5
 * bb = 4 + i8
 * ab = 0.5*(2 + 6) + i0.5*(3 + 7)
 * ba = 0.5*(2 - 6) + i0.5*(3 - 7)
 **/
void splitExchange(int atomNum, int ncao) {
	const std::string cartDict[] = {"x", "y", "z"};

	std::ifstream gxfock("gxfock");
	if (gxfock.fail()) std::cout << "could not read gxfock!" << std::endl;
	std::string line;

	const int caoSize = ncao*(ncao+1)/2;

	for (int c=1; c<=8; c++) {
		getline(gxfock, line);
		getline(gxfock, line);
		getline(gxfock, line);
		for (int I=0; I<atomNum; I++) {
			for (int alpha=0; alpha<=2; alpha++) {
				getline(gxfock, line);
				getline(gxfock, line);
				getline(gxfock, line);

				std::ofstream fil("exch" + std::to_string(I) + cartDict[alpha] + std::to_string(c));
				fil << ncao << "\n";

				for (int i=0; i<caoSize; i++) {
					getline(gxfock, line);
					fil << line << "\n";
				}

				fil.close();
			}
		}
	}
}

void split1efiles(int atomNum, int ncao) {
	/*****************************************************************************
	 *                             overlap and stuff                             *
	 ****************************************************************************/
	std::string line;
	// split sbra and sket files
	std::cout << "\nsplitting sbra and sket files...";
	std::ifstream sbraRe("sbra.r");
	std::ifstream sbraIm("sbra.i");
	std::ifstream sketRe("sket.r");
	std::ifstream sketIm("sket.i");
	std::ifstream jfockRe("gjfock.r");
	std::ifstream jfockIm("gjfock.i");
	if (sbraRe.fail()) std::cout << "\nERROR reading some files!\n";
	if (sbraIm.fail()) std::cout << "\nERROR reading some files!\n";
	if (sketRe.fail()) std::cout << "\nERROR reading some files!\n";
	if (sketIm.fail()) std::cout << "\nERROR reading some files!\n";
	if (jfockRe.fail()) std::cout << "\nERROR reading some files!\n";
	if (jfockIm.fail()) std::cout << "\nERROR reading some files!\n";

	const int caoSize = ncao*(ncao+1)/2;

	getline(sbraRe, line);
	getline(sbraIm, line);
	getline(sketRe, line);
	getline(sketIm, line);
	int size = stoi(line);
	
	for (int i=0; i<atomNum; i++) {
		// x
		std::ofstream brx("b" + std::to_string(i) + "x.r");
		std::ofstream bix("b" + std::to_string(i) + "x.i");
		std::ofstream krx("k" + std::to_string(i) + "x.r");
		std::ofstream kix("k" + std::to_string(i) + "x.i");
		std::ofstream jfxr("jf" + std::to_string(i) + "x.r");
		std::ofstream jfxi("jf" + std::to_string(i) + "x.i");
		brx << size << "\n";
		bix << size << "\n";
		krx << size << "\n";
		kix << size << "\n";
		for (int i=0; i<size*size; i++) {
			getline(sbraRe, line);
			brx << line << "\n";
			getline(sbraIm, line);
			bix << line << "\n";
			getline(sketRe, line);
			krx << line << "\n";
			getline(sketIm, line);
			kix << line << "\n";
		}
		jfxr << ncao << "\n";
		jfxi << ncao << "\n";
		getline(jfockRe, line);
		getline(jfockRe, line);
		getline(jfockRe, line);
		getline(jfockIm, line);
		getline(jfockIm, line);
		getline(jfockIm, line);
		for (int i=0; i<caoSize; i++) {
			getline(jfockRe, line);
			jfxr << line << "\n";
			getline(jfockIm, line);
			jfxi << line << "\n";
		}
		jfxi.close();
		jfxr.close();
		brx.close();
		bix.close();
		krx.close();
		kix.close();

		// y
		std::ofstream bry("b" + std::to_string(i) + "y.r");
		std::ofstream biy("b" + std::to_string(i) + "y.i");
		std::ofstream kry("k" + std::to_string(i) + "y.r");
		std::ofstream kiy("k" + std::to_string(i) + "y.i");
		std::ofstream jfyr("jf" + std::to_string(i) + "y.r");
		std::ofstream jfyi("jf" + std::to_string(i) + "y.i");
		bry << size << "\n";
		biy << size << "\n";
		kry << size << "\n";
		kiy << size << "\n";
		//for (int i=0; i<size; i++) {
		for (int i=0; i<size*size; i++) {
			getline(sbraRe, line);
			bry << line << "\n";
			getline(sbraIm, line);
			biy << line << "\n";
			getline(sketRe, line);
			kry << line << "\n";
			getline(sketIm, line);
			kiy << line << "\n";
		}
		jfyr << ncao << "\n";
		jfyi << ncao << "\n";
		getline(jfockRe, line);
		getline(jfockRe, line);
		getline(jfockRe, line);
		getline(jfockIm, line);
		getline(jfockIm, line);
		getline(jfockIm, line);
		for (int i=0; i<caoSize; i++) {
			getline(jfockRe, line);
			jfyr << line << "\n";
			getline(jfockIm, line);
			jfyi << line << "\n";
		}
		jfyi.close();
		jfyr.close();
		bry.close();
		biy.close();
		kry.close();
		kiy.close();

		// z
		std::ofstream brz("b" + std::to_string(i) + "z.r");
		std::ofstream biz("b" + std::to_string(i) + "z.i");
		std::ofstream krz("k" + std::to_string(i) + "z.r");
		std::ofstream kiz("k" + std::to_string(i) + "z.i");
		std::ofstream jfzr("jf" + std::to_string(i) + "z.r");
		std::ofstream jfzi("jf" + std::to_string(i) + "z.i");
		brz << size << "\n";
		biz << size << "\n";
		krz << size << "\n";
		kiz << size << "\n";
		//for (int i=0; i<size; i++) {
		for (int i=0; i<size*size; i++) {
			getline(sbraRe, line);
			brz << line << "\n";
			getline(sbraIm, line);
			biz << line << "\n";
			getline(sketRe, line);
			krz << line << "\n";
			getline(sketIm, line);
			kiz << line << "\n";
		}
		jfzr << ncao << "\n";
		jfzi << ncao << "\n";
		getline(jfockRe, line);
		getline(jfockRe, line);
		getline(jfockRe, line);
		getline(jfockIm, line);
		getline(jfockIm, line);
		getline(jfockIm, line);
		for (int i=0; i<caoSize; i++) {
			getline(jfockRe, line);
			jfzr << line << "\n";
			getline(jfockIm, line);
			jfzi << line << "\n";
		}
		jfzi.close();
		jfzr.close();
		brz.close();
		biz.close();
		krz.close();
		kiz.close();
	}
	jfockRe.close();
	jfockIm.close();
	sbraRe.close();
	sbraIm.close();
	sketRe.close();
	sketIm.close();

	std::cout << "   done\n";


	/*****************************************************************************
	 *                         core hamilton and stuff                           *
	 ****************************************************************************/
	// split hgrad files
	std::cout << "\nsplitting hgrad files...";
	std::ifstream hgradRe("hgrad.r");
	std::ifstream hgradIm("hgrad.i");
	if (hgradRe.fail()) std::cout << "\nERROR reading some files!\n";
	if (hgradIm.fail()) std::cout << "\nERROR reading some files!\n";

	getline(hgradRe, line);
	getline(hgradIm, line);
	size = stoi(line);
	
	for (int i=0; i<atomNum; i++) {
		// x
		std::ofstream hrx("h" + std::to_string(i) + "x.r");
		std::ofstream hix("h" + std::to_string(i) + "x.i");
		hrx << size << "\n";
		hix << size << "\n";
		for (int i=0; i<size; i++) {
			getline(hgradRe, line);
			hrx << line << "\n";
			getline(hgradIm, line);
			hix << line << "\n";
		}
		hrx.close();
		hix.close();

		// y
		std::ofstream hry("h" + std::to_string(i) + "y.r");
		std::ofstream hiy("h" + std::to_string(i) + "y.i");
		hry << size << "\n";
		hiy << size << "\n";
		for (int i=0; i<size; i++) {
			getline(hgradRe, line);
			hry << line << "\n";
			getline(hgradIm, line);
			hiy << line << "\n";
		}
		hry.close();
		hiy.close();

		// z
		std::ofstream hrz("h" + std::to_string(i) + "z.r");
		std::ofstream hiz("h" + std::to_string(i) + "z.i");
		hrz << size << "\n";
		hiz << size << "\n";
		for (int i=0; i<size; i++) {
			getline(hgradRe, line);
			hrz << line << "\n";
			getline(hgradIm, line);
			hiz << line << "\n";
		}
		hrz.close();
		hiz.close();
	}
	hgradRe.close();
	hgradIm.close();

	std::cout << "   done\n\n\n";
}

Eigen::VectorXcd berryRHS(const int nuc, const int cart, const Eigen::VectorXcd& ailkasym) {
	//std::cout << "BERRY\n\n" << std::flush;
	std::cout << std::setprecision(10);
	const std::string cartDict[] = {"x", "y", "z"};

	//const int nuc = std::stoi(argv[1]);
	//const int cart = std::stoi(argv[2]);

	std::cout << " :: calculating RHS for atom " << nuc+1 << " " << cartDict[cart] << "\n" << std::flush;
	//const int nuc = 0;
	//const int cart = 2; // 0=x  1=y  2=z
	//std::cout << "of nuclues:   " << nuc << "        (starting at 0)\n";
	//std::cout << "cart = " << cart << "    (0=x  1=y  2=z)\n";
	//std::cout << "\n\n";
	

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
	std::vector<int> g(atomNum);
	int sSize = 0;
	int pSize = 0;
	int dSize = 0;
	int fSize = 0;
	int gSize = 0;
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
							for (; i<word.length()-1; i++) { // f
								if (word[i] == 'f') break;
								f[atom] *= 10;
								f[atom] += word[i] - '0';
							}
							fSize += f[atom];
							i++;
							g.push_back(0);
							for (; i<word.length()-1; i++) { // g
								if (word[i] == 'f') break;
								g[atom] *= 10;
								g[atom] += word[i] - '0';
							}
							gSize += g[atom];
						}
					}
				}
			}
		}
	}
	std::vector<int> basisFunctionLength(atomNum);
	int matrixSize = 0;
	for (int i=0; i<atomNum; i++) {
		basisFunctionLength[i] = s[i] + 3*p[i] + 5*d[i] + 7*f[i] + 9*g[i];
		matrixSize += s[i] + 3*p[i] + 5*d[i] + 7*f[i] + 9*g[i];
	}



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


	// read smat just for fun and spin zeeman...
	//std::cout << "reading smat...";
	const auto smat = readHerm("smat");
	//std::cout << "   done\n";
	//std::cout << "\nsmat:\n" << smat << "\n\n";
	
	/**********************************************************************
	 *                         read control file                          *
	 **********************************************************************
	 * 
	 * current infos taken from control file:
	 * 	- number of occupied spinors
	 */

	// looking for occupied spinors in control file
	std::cout << "\treading control file for some additional infos..." << std::flush;
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
	if (Bnorm == 0) { Bx=0; By=0; Bz=0; }
	// TESLA TODO ÄNDERN HIER ACHTUNG, PLS DONT FORGET TO COMMENT THIS OUT, FUTURE IDIOT <=====================================================================================================
	//Bx /= 2.35051756758e5;
	//By /= 2.35051756758e5;
	//Bz /= 2.35051756758e5;
	//std::cout << "scaled B vector: " << Bx << "  " << By << "  " << Bz << "\n";
	//std::cout << "Bnorm = " << Bnorm << "\n";

	// end of control file section
	control.close();
	std::cout << "   done\n" << std::flush;

	
	std::cout << "\treading 1e matrices..." << std::flush;
	const auto bra = readMatrixTransform("b" + std::to_string(nuc) + cartDict[cart]);
	const auto ket = readMatrixTransform("k" + std::to_string(nuc) + cartDict[cart]);
	const auto snx = (bra + ket).conjugate();
	const auto hnx = readHerm("h" + std::to_string(nuc) + cartDict[cart]);
	std::cout << "   done\n" << std::flush;


	/*****************************************************************************
	 *                             calculate fock matrix                         *
	 ****************************************************************************/
	std::cout << "\tcalculating density matrix..." << std::flush;
	std::vector<double> epsilon;
	auto spinor = readSpinor(epsilon);
	int spinorSize = spinor.rows();
	matrixSize = spinorSize/2;
	Eigen::MatrixXcd denMat = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	for (int k=0; k<spinorSize; k++) {
		for (int l=0; l<spinorSize; l++) {
			for (int i=0; i<nocc; i++) {
				denMat(k, l) += std::conj(spinor(k, i)) * spinor(l, i);
			}
		}
	}
	std::cout << "   done\n" << std::flush;


	// ==================================== SPIN ZEEMAN ===================================

	// spin-Zeeman contribution (SAO)
	std::cout << "\tcalculating spin Zeeman contribution" << std::flush;
	Eigen::MatrixXcd zeemanx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	Eigen::MatrixXcd zeemany = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	Eigen::MatrixXcd zeemanz = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	//fSAOZeeman << 0.5*Bz*smat, 0.5*std::complex<double>(Bx, -By)*smat,
	//	      0.5*std::complex<double>(Bx, By)*smat, -0.5*Bz*smat;
	zeemanx << Eigen::MatrixXcd::Zero(matrixSize, matrixSize), 0.5*Bx*smat,
		0.5*Bx*smat, Eigen::MatrixXcd::Zero(matrixSize, matrixSize);
	zeemany << Eigen::MatrixXcd::Zero(matrixSize, matrixSize), std::complex<double>(0, 0.5)*By*smat,
		std::complex<double>(0, -0.5)*By*smat, Eigen::MatrixXcd::Zero(matrixSize, matrixSize);
	zeemanz << 0.5*Bz*smat, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), -0.5*Bz*smat;
	auto fockSAO = zeemanx;// + zeemany + zeemanz;
	fockSAO += zeemany;
	fockSAO += zeemanz;
	std::cout << "   done\n" << std::flush;

	std::cout << "\tcalculating derivative of spin Zeeman contribution" << std::flush;
	constexpr double spinFactor = 1.0;

	// transform to spinor basis
	zeemanx = spinFactor * spinor.adjoint() * zeemanx * spinor;
	zeemany = spinFactor * spinor.adjoint() * zeemany * spinor;
	zeemanz = spinFactor * spinor.adjoint() * zeemanz * spinor;
	//std::cout << "F Zeeman:\n" << fZeeman << "\n\n";
	std::complex<double> eSpinZeemanX = 0.0;
	std::complex<double> eSpinZeemanY = 0.0;
	std::complex<double> eSpinZeemanZ = 0.0;
	for (int i=0; i<nocc; i++) {
		eSpinZeemanX += zeemanx(i, i);
		eSpinZeemanY += zeemany(i, i);
		eSpinZeemanZ += zeemanz(i, i);
	}


	// derivative of spin-Zeeman contribution (SAO)
	Eigen::MatrixXcd zxnx(spinorSize, spinorSize);
	Eigen::MatrixXcd zynx(spinorSize, spinorSize);
	Eigen::MatrixXcd zznx(spinorSize, spinorSize);
	zxnx << Eigen::MatrixXcd::Zero(matrixSize, matrixSize), 0.5*Bx*snx,
		0.5*Bx*snx, Eigen::MatrixXcd::Zero(matrixSize, matrixSize);
	zynx << Eigen::MatrixXcd::Zero(matrixSize, matrixSize), std::complex<double>(0, -0.5)*By*snx,
		std::complex<double>(0, +0.5)*By*snx, Eigen::MatrixXcd::Zero(matrixSize, matrixSize);
	zznx << 0.5*Bz*snx, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), -0.5*Bz*snx;
	zxnx *= spinFactor;
	zynx *= spinFactor;
	zznx *= spinFactor;

	//std::cout << std::setprecision(3);
	//std::cout << "\n\nzz nx:\n" << zznx << "\n";

	Eigen::MatrixXcd fnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	fnx =  zxnx.transpose();
	fnx += zynx.transpose();
	fnx += zznx.transpose();
	//std::cout << "\nFZ nx total:\n" << fnx << "\n\n";
	Eigen::MatrixXcd zetotnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	zetotnx = zxnx;
	zetotnx += zynx;
	zetotnx += zznx;
	std::cout << "   done\n" << std::flush;


	// fourcenter derivatives in SAO basis
	const auto beginread = std::chrono::high_resolution_clock::now();
	std::cout << std::flush << "\treading g4c..." << std::flush;
	//const auto fcinx = fciAlt(nuc, cart);
	//const auto endread = std::chrono::high_resolution_clock::now();
	//const auto elapsedread = std::chrono::duration_cast<std::chrono::milliseconds>(endread-beginread);
	//std::cout << "\tdone reading after " << elapsedread.count()*1e-3 << " s\n\n";
	//std::cout << "\tcalculating 2e part of gradient...\n" << std::flush;

	//auto begin1 = std::chrono::high_resolution_clock::now();
	//std::cout << "\tdefine spin components... " << std::flush;
	//// double den kack
	//fourD fcinxD(spinorSize,
	//		std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
	//			std::vector<std::vector<std::complex<double>>>(spinorSize,
	//				std::vector<std::complex<double>>(spinorSize))));

	//for (int i=0; i<spinorSize; i++) {
	//	for (int j=0; j<spinorSize; j++) {
	//		for (int k=0; k<spinorSize; k++) {
	//			for (int l=0; l<spinorSize; l++) {
	//				fcinxD[i][j][k][l] = (0, 0);
	//			}
	//		}
	//	}
	//}
	//for (int i=0; i<spinorSize/2; i++) {
	//	for (int j=0; j<spinorSize/2; j++) {
	//		for (int k=0; k<spinorSize/2; k++) {
	//			for (int l=0; l<spinorSize/2; l++) {
	//				fcinxD[i][j][k][l] = fcinx[i][j][k][l];
	//				fcinxD[i][j][k+spinorSize/2][l+spinorSize/2] = fcinx[i][j][k][l];
	//				fcinxD[i+spinorSize/2][j+spinorSize/2][k][l] = fcinx[i][j][k][l];
	//				fcinxD[i+spinorSize/2][j+spinorSize/2][k+spinorSize/2][l+spinorSize/2] = fcinx[i][j][k][l];
	//			}
	//		}
	//	}
	//}
	//auto end1 = std::chrono::high_resolution_clock::now();
	//auto elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1-begin1);
	//printf(" done after %.3fs\n", elapsed1.count()*1e-3);

	

	// ==================================== core hamilton ===================================

	// direkt in spinor basis
	Eigen::MatrixXcd hmatBig(spinorSize, spinorSize);
	const auto hmat = readHerm("hmat");
	hmatBig << hmat, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), hmat;
	fockSAO += hmatBig;
	const auto hmatkek = hmatBig;
	hmatBig = spinor.adjoint() * hmatBig * spinor;

	// derivative
	Eigen::MatrixXcd hnxBig(spinorSize, spinorSize);
	hnxBig << hnx.conjugate(), Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), hnx.conjugate();
	fnx += hnxBig.transpose();




	// ==================================== two electron matrix ===================================
	//
	// calculate Fock matrix in spinor basis
	// F = h + G + ZF
	// calc G first
	Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	Eigen::MatrixXcd K = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	//Eigen::MatrixXcd Cnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	//Eigen::MatrixXcd Knx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	auto begin2 = std::chrono::high_resolution_clock::now();
	std::cout << "\tcalculating derivative of G..." << std::flush;
	std::cout << "\n";
	//for (int k=0; k<spinorSize; k++) {
	//	for (int l=0; l<spinorSize; l++) {
	//		for (int m=0; m<spinorSize; m++) {
	//			for (int n=0; n<spinorSize; n++) {
	//				Cnx(k, l) += denMat(n, m) * fcinxD[k][l][m][n];
	//				Knx(k, l) -= denMat(n, m) * fcinxD[k][n][m][l];
	//			}
	//		}
	//	}
	//}
	//std::cout << "\n";

	// read Coulomb derivative
	const auto Cnx = readCoulomb("jf" + std::to_string(nuc) + cartDict[cart]);
	//auto diff = Jnxread - Cnx.block(0, 0, spinorSize/2, spinorSize/2);
	//double diffmax = -1.0;
	//for (int i=0; i<spinorSize/2; i++) {
	//	for (int j=0; j<spinorSize/2; j++) {
	//		diffmax = std::max(diffmax, abs(diff(i, j)));
	//	}
	//}
	//std::cout << " --------------->  diffmax coul = " << diffmax << "\n";

	//std::cout << "\ndJ/dIa real:\n" << std::fixed << std::setprecision(5) << Cnx.real() << "\n";
	//std::cout << "\nJnx read real:\n" << std::fixed << std::setprecision(5) << Jnxread.real() << "\n";
	//std::cout << "\ndJ/dIa imag:\n" << std::fixed << std::setprecision(5) << Cnx.imag() << "\n";
	//std::cout << "\nJnx read imag:\n" << std::fixed << std::setprecision(5) << Jnxread.imag() << "\n";

	const auto Knx1 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 1);
	const auto Knx2 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 2);
	const auto Knx3 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 3);
	const auto Knx4 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 4);
	const auto Knx5 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 5);
	const auto Knx6 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 6);
	const auto Knx7 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 7);
	const auto Knx8 = readExchange("exch" + std::to_string(nuc) + cartDict[cart], 8);

	Eigen::MatrixXcd Knx(spinorSize, spinorSize);
	Knx.block(0, 0, spinorSize/2, spinorSize/2).real() = Knx1;				// Re aa
	Knx.block(0, 0, spinorSize/2, spinorSize/2).imag() = Knx5;				// Im aa
	Knx.block(spinorSize/2, spinorSize/2, spinorSize/2, spinorSize/2).real() = Knx4;	// Re bb
	Knx.block(spinorSize/2, spinorSize/2, spinorSize/2, spinorSize/2).imag() = Knx8;	// Im bb
	Knx.block(0, spinorSize/2, spinorSize/2, spinorSize/2).real() = 0.5*(Knx2 - Knx6);	// Re ab
	Knx.block(0, spinorSize/2, spinorSize/2, spinorSize/2).imag() = -0.5*(Knx3 - Knx7);	// Im ab
	Knx.block(spinorSize/2, 0, spinorSize/2, spinorSize/2).real() = 0.5*(Knx2 + Knx6);	// Re ba
	Knx.block(spinorSize/2, 0, spinorSize/2, spinorSize/2).imag() = 0.5*(Knx3 + Knx7);	// Im ba
												
	//std::cout << "\ndK/dIa real:\n" << std::fixed << std::setprecision(5) << Knx.real() << "\n";
	//std::cout << "\ndK/dIa real read:\n" << std::fixed << std::setprecision(5) << Knxread.real() << "\n";
	//std::cout << "\ndK/dIa imag:\n" << std::fixed << std::setprecision(5) << Knx.imag() << "\n";
	//std::cout << "\ndK/dIa imag read:\n" << std::fixed << std::setprecision(5) << Knxread.imag() << "\n";
	
	//const auto diffexaa = (Knx - Knxread).block(0, 0, spinorSize/2, spinorSize/2);
	//const auto diffexbb = (Knx - Knxread).block(spinorSize/2, spinorSize/2, spinorSize/2, spinorSize/2);
	//const auto diffexab = (Knx - Knxread).block(0, spinorSize/2, spinorSize/2, spinorSize/2);
	//const auto diffexba = (Knx - Knxread).block(spinorSize/2, 0, spinorSize/2, spinorSize/2);
	//double diffmaxaa = -1.0;
	//double diffmaxbb = -1.0;
	//double diffmaxab = -1.0;
	//double diffmaxba = -1.0;
	//for (int i=0; i<spinorSize/2; i++) {
	//	for (int j=0; j<spinorSize/2; j++) {
	//		diffmaxaa = std::max(diffmaxaa, abs(diffexaa(i, j)));
	//		diffmaxbb = std::max(diffmaxbb, abs(diffexbb(i, j)));
	//		diffmaxab = std::max(diffmaxab, abs(diffexab(i, j)));
	//		diffmaxba = std::max(diffmaxba, abs(diffexba(i, j)));
	//	}
	//}
	//std::cout << " --------------->  diffmax exch aa = " << diffmaxaa << "\n";
	//std::cout << " --------------->  diffmax exch bb = " << diffmaxbb << "\n";
	//std::cout << " --------------->  diffmax exch ab = " << diffmaxab << "\n";
	//std::cout << " --------------->  diffmax exch ba = " << diffmaxba << "\n";

	//std::cout << std::fixed << std::setprecision(8) << "Knx real:\n" << Knx.block(0, spinorSize/2, spinorSize/2, spinorSize/2).real() << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "Knx read real:\n" << Knxread.block(0, spinorSize/2, spinorSize/2, spinorSize/2).real() << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "Knx imag:\n" << Knx.block(0, spinorSize/2, spinorSize/2, spinorSize/2).imag() << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "Knx read imag:\n" << Knxread.block(0, spinorSize/2, spinorSize/2, spinorSize/2).imag() << "\n\n";
	std::cout << std::defaultfloat;


	auto end2 = std::chrono::high_resolution_clock::now();
	auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2-begin2);
	printf(" done after %.3fs\n", elapsed2.count()*1e-3);
	std::cout << std::flush;

	//std::cout << "density real:\n" << denMat.real() << "\n\n";
	//std::cout << "density imag:\n" << denMat.imag() << "\n\n";
	//std::cout << "coulomb:\n" << C << "\n\n";
	//std::cout << "exchange:\n" << K << "\n\n";
	//fockSAO += C.transpose();
	//fockSAO += K.transpose();
	std::cout << "\tupdating fock derivative..." << std::flush;
	fnx.block(0, 0, spinorSize/2, spinorSize/2) += Cnx.transpose();
	fnx.block(spinorSize/2, spinorSize/2, spinorSize/2, spinorSize/2) += Cnx.transpose();
	fnx += Knx.transpose();
	std::cout << " done.\n" << std::flush;



	/*****************************************************************************
	 *                              gradient stuff                               *
	 ****************************************************************************/
	//std::cout << "calculating gradient...\n";
	Eigen::MatrixXcd snxBig(spinorSize, spinorSize);
	snxBig << snx, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), snx;

	// transpose fnx
	fnx.transposeInPlace();


	const size_t nvirt = spinorSize - nocc;
	Eigen::VectorXcd b0ai = Eigen::VectorXcd::Zero(nocc*nvirt);

	const auto fnx2cMO = spinor.adjoint() * fnx.conjugate() * spinor;
	const auto snx2cMO = spinor.adjoint() * snxBig.conjugate() * spinor;


	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nfnx AO:\n" << fnx << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nfnx MO:\n" << fnx2cMO << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nsnxAO:\n" << snx << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nsnx2cMO real:\n" << snx2cMO.real() << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nsnx2cMO imag:\n" << snx2cMO.imag() << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nc alter:\n" << Calternative << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nKspinor:\n" << Kspinor << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nK alter:\n" << Kalternative << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nhmat + c + k spinor:\n" << hmatBig.transpose() + Calternative + Kalternative << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nCoul AO:\n" << spinor.adjoint().inverse() * Calternative.transpose() * spinor.inverse() << "\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nCoul AO:\n" << C << "\n";
	//std::cout << "\n\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nExch AO:\n" << spinor.adjoint().inverse() * Kalternative.transpose() * spinor.inverse() << "\n";
	//std::cout << std::fixed << std::setprecision(8) << "\nExch AO:\n" << K << "\n";
	//std::cout << std::defaultfloat;



	// constructing snxVector
	Eigen::VectorXcd snxVec(nocc*nocc);
	for (int k=0; k<nocc; k++) {
		for (int l=0; l<nocc; l++) {
			snxVec(l + nocc*k) = snx2cMO(k, l);
		}
	}

	auto begin = std::chrono::high_resolution_clock::now();
	std::cout << "\tputting RHS together..." << std::flush;
	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			const int index = i*nvirt + a - nocc;
			//b0ai(index) = -fnx2cMO(a, i);
			//b0ai(index) += snx2cMO(a, i) * epsilon[i];
		}
	}
	std::cout << " 1e part done" << std::flush;
	for (int i=0; i<nocc; i++) {
		for (int k=0; k<nocc; k++) {
			for (int l=0; l<nocc; l++) {
				for (int a=nocc; a<spinorSize; a++) {
					b0ai(i*nvirt + a - nocc) += snxVec(l+nocc*k) * ailkasym( (a-nocc) + nvirt*l + nvirt*nocc*k + nvirt*nocc*nocc*i );
				}
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	printf("\ttotal done after %.3fs\n", elapsed.count()*1e-3);
	std::cout << std::flush;

	// write b0ai to file
	std::ofstream outRe("b0ai" + std::to_string(nuc) + "_" + std::to_string(cart) + ".r");
	std::ofstream outIm("b0ai" + std::to_string(nuc) + "_" + std::to_string(cart) + ".i");
	outRe << nocc*nvirt << "\n";
	outIm << nocc*nvirt << "\n";
	for (int i=0; i<nocc*nvirt; i++) {
		outRe << std::setprecision(12) <<  b0ai(i).real() << "\n";
		outIm << std::setprecision(12) <<  b0ai(i).imag() << "\n";
	}
	outRe.close();
	outIm.close();

	//std::cout << "b0ai:\n" << b0ai << "\n";


	std::cout << "done calculating RHS\n" << std::flush;
	return b0ai;
}
