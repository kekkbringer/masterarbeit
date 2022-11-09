#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <complex>
#include <fstream>
#include <iomanip>
#include <string>

#include "misc.hpp"
#include "cphf.hpp"
#include "read_spinor.hpp"
#include "read_herm.hpp"
#include "read_fourcenter.hpp"
#include "calc_stabmat.hpp"
#include "fci_grad.hpp"

#define DEBUG 1
#define IFDBG if constexpr (DEBUG)

#define MAJOR_VERSION 0
#define MINOR_VERSION 0
#define PATCH_VERSION 1

int main(int argc, char* argv[]) {
	std::cout << "BERRY\n\n" << std::flush;
	std::cout << std::setprecision(10);

	const int nuc = std::stoi(argv[1]);
	const int cart = std::stoi(argv[2]);

	std::cout << "currently im only calculating the gradient (dE/dN_x)\n";
	//const int nuc = 0;
	//const int cart = 2; // 0=x  1=y  2=z
	std::cout << "of nuclues:   " << nuc << "        (starting at 0)\n";
	std::cout << "cart = " << cart << "    (0=x  1=y  2=z)\n";
	std::cout << "\n\n";
	

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
	std::cout << "matrix size: " << matrixSize << "\n\n";



	std::cout << "      number of atoms: " << atomNum << "\n";
	int lower = 0;
	int upper = 0;
	for (int i=0; i<atomNum; i++)  {
		std::cout << atoms[i] << "   " << basisFunctionLength[i] << "\n";
		if (i<nuc) lower += basisFunctionLength[i];
		if (i<=nuc) upper += basisFunctionLength[i];
	}
	std::cout << "lower = " << lower << "\n";
	std::cout << "upper = " << upper << "\n";


	// read smat just for fun and spin zeeman...
	std::cout << "reading smat...\n";
	const auto smat = readHerm("smat");
	//std::cout << "\nsmat:\n" << smat << "\n\n";
	
	/**********************************************************************
	 *                         read control file                          *
	 **********************************************************************
	 * 
	 * current infos taken from control file:
	 * 	- number of occupied spinors
	 */

	// looking for occupied spinors in control file
	std::cout << "\n:: reading control file for some additional infos...\n" << std::flush;
	int nocc = 0;
	double Bx = 0.0;
	double By = 0.0;
	double Bz = 0.0;
	double Bnorm = 0.0;
	std::ifstream control("control");
	if (control.fail()) {std::cout << "\nWARNING: cant read control file!\n"; return -1;}
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
	std::cout << "first " << nocc << " spinors are occupied\n";
	std::cout << "magnetic field strength: " << Bnorm << "\n";
	std::cout << "  direction:             " << Bx << "  " << By << "  " << Bz << "\n";
	// normalize B
	const double lambda = Bnorm / sqrt(Bx*Bx + By*By + Bz*Bz);
	Bx *= lambda;
	By *= lambda;
	Bz *= lambda;
	std::cout << "scaled B vector: " << Bx << "  " << By << "  " << Bz << "\n";

	// end of control file section
	control.close();

	std::cout << "\n" << std::flush;



	/*****************************************************************************
	 *                             overlap and stuff                             *
	 ****************************************************************************/
	// split sbra and sket files
	std::cout << "\nsplitting sbra and sket files...\n";
	std::ifstream sbraRe("sbra.r");
	std::ifstream sbraIm("sbra.i");
	std::ifstream sketRe("sket.r");
	std::ifstream sketIm("sket.i");
	if (sbraRe.fail()) std::cout << "\nERROR reading some files!\n";
	if (sbraIm.fail()) std::cout << "\nERROR reading some files!\n";
	if (sketRe.fail()) std::cout << "\nERROR reading some files!\n";
	if (sketIm.fail()) std::cout << "\nERROR reading some files!\n";

	getline(sbraRe, line);
	getline(sbraIm, line);
	getline(sketRe, line);
	getline(sketIm, line);
	const int size = stoi(line);
	
	for (int i=0; i<atomNum; i++) {
		// x
		std::ofstream brx("b" + std::to_string(i) + "x.r");
		std::ofstream bix("b" + std::to_string(i) + "x.i");
		std::ofstream krx("k" + std::to_string(i) + "x.r");
		std::ofstream kix("k" + std::to_string(i) + "x.i");
		brx << size << "\n";
		bix << size << "\n";
		krx << size << "\n";
		kix << size << "\n";
		for (int i=0; i<size; i++) {
			getline(sbraRe, line);
			brx << line << "\n";
			getline(sbraIm, line);
			bix << line << "\n";
			getline(sketRe, line);
			krx << line << "\n";
			getline(sketIm, line);
			kix << line << "\n";
		}
		brx.close();
		bix.close();
		krx.close();
		kix.close();

		// y
		std::ofstream bry("b" + std::to_string(i) + "y.r");
		std::ofstream biy("b" + std::to_string(i) + "y.i");
		std::ofstream kry("k" + std::to_string(i) + "y.r");
		std::ofstream kiy("k" + std::to_string(i) + "y.i");
		bry << size << "\n";
		biy << size << "\n";
		kry << size << "\n";
		kiy << size << "\n";
		for (int i=0; i<size; i++) {
			getline(sbraRe, line);
			bry << line << "\n";
			getline(sbraIm, line);
			biy << line << "\n";
			getline(sketRe, line);
			kry << line << "\n";
			getline(sketIm, line);
			kiy << line << "\n";
		}
		bry.close();
		biy.close();
		kry.close();
		kiy.close();

		// z
		std::ofstream brz("b" + std::to_string(i) + "z.r");
		std::ofstream biz("b" + std::to_string(i) + "z.i");
		std::ofstream krz("k" + std::to_string(i) + "z.r");
		std::ofstream kiz("k" + std::to_string(i) + "z.i");
		brz << size << "\n";
		biz << size << "\n";
		krz << size << "\n";
		kiz << size << "\n";
		for (int i=0; i<size; i++) {
			getline(sbraRe, line);
			brz << line << "\n";
			getline(sbraIm, line);
			biz << line << "\n";
			getline(sketRe, line);
			krz << line << "\n";
			getline(sketIm, line);
			kiz << line << "\n";
		}
		brz.close();
		biz.close();
		krz.close();
		kiz.close();
	}
	sbraRe.close();
	sbraIm.close();
	sketRe.close();
	sketIm.close();


	// construct S_ab^(N_x)
	// N is in nuc
	// x is in cart
	const std::string cartDict[] = {"x", "y", "z"};
	const auto bra = readHerm("b" + std::to_string(nuc) + cartDict[cart]);
	const auto ket = readHerm("k" + std::to_string(nuc) + cartDict[cart]);
	const auto snx = (bra + ket).conjugate();
	//std::cout << "snx dim: " << snx.rows() << " x " << snx.cols() << "\n";
	//std::cout << "\nsnx:\n" << snx << "\n\n";
	//*/

	/*****************************************************************************
	 *                         core hamilton and stuff                           *
	 ****************************************************************************/
	// split hgrad files
	std::cout << "\nsplitting hgrad files...\n";
	std::ifstream hgradRe("hgrad.r");
	std::ifstream hgradIm("hgrad.i");
	if (hgradRe.fail()) std::cout << "\nERROR reading some files!\n";
	if (hgradIm.fail()) std::cout << "\nERROR reading some files!\n";

	getline(hgradRe, line);
	getline(hgradIm, line);
	//const int size = stoi(line);
	
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


	// construct H_ab^(N_x)
	// N is in nuc
	// x is in cart
	const auto hnx = readHerm("h" + std::to_string(nuc) + cartDict[cart]);
	//std::cout << "hnx dim: " << hnx.rows() << " x " << hnx.cols() << "\n";
	//std::cout << "\nhnx:\n" << hnx << "\n\n";
	//*/

	/*****************************************************************************
	 *                             calculate fock matrix                         *
	 ****************************************************************************/
	std::cout << "calculating fock matrix...\n";
	std::vector<double> epsilon;
	auto spinor = readSpinor(epsilon);
	int spinorSize = spinor.rows();
	Eigen::MatrixXcd denMat = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	for (int k=0; k<spinorSize; k++) {
		for (int l=0; l<spinorSize; l++) {
			for (int i=0; i<nocc; i++) {
				denMat(k, l) += std::conj(spinor(k, i)) * spinor(l, i);
			}
		}
	}


	// ==================================== SPIN ZEEMAN ===================================

	// spin-Zeeman contribution (SAO)
	std::cout << "calculating spin Zeeman contribution\n";
	Eigen::MatrixXcd zeemanx(spinorSize, spinorSize);
	Eigen::MatrixXcd zeemany(spinorSize, spinorSize);
	Eigen::MatrixXcd zeemanz(spinorSize, spinorSize);
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

	std::cout << "calculating derivative of spin Zeeman contribution\n";
	const double spinFactor = 1.0;

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
	zynx << Eigen::MatrixXcd::Zero(matrixSize, matrixSize), std::complex<double>(0, 0.5)*By*snx,
		std::complex<double>(0, -0.5)*By*snx, Eigen::MatrixXcd::Zero(matrixSize, matrixSize);
	zznx << 0.5*Bz*snx, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), -0.5*Bz*snx;
	zxnx *= spinFactor;
	zynx *= -1.0 * spinFactor;
	zznx *= spinFactor;

	//std::cout << std::setprecision(3);
	//std::cout << "\n\nsnx:\n" << snx << "\n";

	Eigen::MatrixXcd fnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	fnx =  zxnx.transpose();
	fnx += zynx.transpose();
	fnx += zznx.transpose();
	Eigen::MatrixXcd zetotnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	zetotnx = zxnx;
	zetotnx += zynx;
	zetotnx += zznx;

	//std::cout << "\n\ndFock-Zeeman-SAO-x/dNx:\n" << zxnx << "\n";

	

	// ==================================== fourcenter ===================================

	// double dim of fci
	auto fci = readFourcenter();
	//auto fci = fciAlt(nuc, cart, lower, upper);
	fourD fciD(spinorSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
				std::vector<std::vector<std::complex<double>>>(spinorSize,
					std::vector<std::complex<double>>(spinorSize))));

	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			for (int k=0; k<spinorSize; k++) {
				for (int l=0; l<spinorSize; l++) {
					fciD[i][j][k][l] = (0, 0);
				}
			}
		}
	}
	for (int i=0; i<spinorSize/2; i++) {
		for (int j=0; j<spinorSize/2; j++) {
			for (int k=0; k<spinorSize/2; k++) {
				for (int l=0; l<spinorSize/2; l++) {
					fciD[i][j][k][l] = fci[i][j][k][l];
					fciD[i][j][k+spinorSize/2][l+spinorSize/2] = fci[i][j][k][l];
					fciD[i+spinorSize/2][j+spinorSize/2][k][l] = fci[i][j][k][l];
					fciD[i+spinorSize/2][j+spinorSize/2][k+spinorSize/2][l+spinorSize/2] = fci[i][j][k][l];
				}
			}
		}
	}



	// fourcenter derivatives in SAO basis
	std::cout << std::flush << "calculating/reading g4c...\n" << std::flush;
	const auto fcinx = fciAlt(nuc, cart, lower, upper);
	std::cout << "done reading!\n\ncalculating 2e part of gradient...\n" << std::flush;

	// double den kack
	fourD fcinxD(spinorSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
				std::vector<std::vector<std::complex<double>>>(spinorSize,
					std::vector<std::complex<double>>(spinorSize))));

	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			for (int k=0; k<spinorSize; k++) {
				for (int l=0; l<spinorSize; l++) {
					fcinxD[i][j][k][l] = (0, 0);
				}
			}
		}
	}
	for (int i=0; i<spinorSize/2; i++) {
		for (int j=0; j<spinorSize/2; j++) {
			for (int k=0; k<spinorSize/2; k++) {
				for (int l=0; l<spinorSize/2; l++) {
					fcinxD[i][j][k][l] = fcinx[i][j][k][l];
					fcinxD[i][j][k+spinorSize/2][l+spinorSize/2] = fcinx[i][j][k][l];
					fcinxD[i+spinorSize/2][j+spinorSize/2][k][l] = fcinx[i][j][k][l];
					fcinxD[i+spinorSize/2][j+spinorSize/2][k+spinorSize/2][l+spinorSize/2] = fcinx[i][j][k][l];
				}
			}
		}
	}


	

	// ==================================== core hamilton ===================================

	// direkt in spinor basis
	Eigen::MatrixXcd hmatBig(spinorSize, spinorSize);
	const auto hmat = readHerm("hmat");
	hmatBig << hmat, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), hmat;
	fockSAO += hmatBig;
	const auto hmatkek = hmatBig;
	hmatBig = spinor.adjoint() * hmatBig * spinor;

	// calculate 1e energy
	std::complex<double> e1 = 0;
	for (int a=0; a<nocc; a++) {
		e1 += hmatBig(a, a);
	}
	

	// derivative
	Eigen::MatrixXcd hnxBig(spinorSize, spinorSize);
	hnxBig << hnx.conjugate(), Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), hnx.conjugate();
	fnx += hnxBig.transpose();





	// ==================================== nuc-nuc interaction ===================================

	// calculate electric potential energy of nuclei
	double enuc = 0.0;
	for (int i=0; i<atomNum; i++) {
		for (int j=0; j<atomNum; j++) {
			if (i==j) continue;
			enuc += chargeOf(atoms[i]) * chargeOf(atoms[j])
				/ sqrt(  (coordx[i]-coordx[j])*(coordx[i]-coordx[j])
				   + (coordy[i]-coordy[j])*(coordy[i]-coordy[j])
				   + (coordz[i]-coordz[j])*(coordz[i]-coordz[j])  );
		}
	}
	enuc *= 0.5;


	// derivatives
	double enucnx[atomNum][3];
	// x derivatives
	for (int a=0; a<atomNum; a++) {
		enucnx[a][0] = 0.0;
		for (int b=0; b<atomNum; b++) {
			if (a==b) continue;
			double rab3 = pow( (coordx[a]-coordx[b])*(coordx[a]-coordx[b])
				   + (coordy[a]-coordy[b])*(coordy[a]-coordy[b])
				   + (coordz[a]-coordz[b])*(coordz[a]-coordz[b]), 1.5);
			enucnx[a][0] += chargeOf(atoms[b]) * (coordx[b] - coordx[a]) / rab3;
		}
		enucnx[a][0] *= chargeOf(atoms[a]);
	}
	// y derivatives
	for (int a=0; a<atomNum; a++) {
		enucnx[a][1] = 0.0;
		for (int b=0; b<atomNum; b++) {
			if (a==b) continue;
			double rab3 = pow( (coordx[a]-coordx[b])*(coordx[a]-coordx[b])
				   + (coordy[a]-coordy[b])*(coordy[a]-coordy[b])
				   + (coordz[a]-coordz[b])*(coordz[a]-coordz[b]), 1.5);
			enucnx[a][1] += chargeOf(atoms[b]) * (coordy[b] - coordy[a]) / rab3;
		}
		enucnx[a][1] *= chargeOf(atoms[a]);
	}
	// z derivatives
	for (int a=0; a<atomNum; a++) {
		enucnx[a][2] = 0.0;
		for (int b=0; b<atomNum; b++) {
			if (a==b) continue;
			double rab3 = pow( (coordx[a]-coordx[b])*(coordx[a]-coordx[b])
				   + (coordy[a]-coordy[b])*(coordy[a]-coordy[b])
				   + (coordz[a]-coordz[b])*(coordz[a]-coordz[b]), 1.5);
			enucnx[a][2] += chargeOf(atoms[b]) * (coordz[b] - coordz[a]) / rab3;
		}
		enucnx[a][2] *= chargeOf(atoms[a]);
	}

	// print derivatives d(Vnn)/dNx
	std::cout << "\nderivatives of Vnn:\n";
	for (int a=0; a<atomNum; a++) {
		std::cout << enucnx[a][0] << "   " << enucnx[a][1] << "   " << enucnx[a][2] << "\n";
	}




	// ==================================== small debug print ===================================

	std::cout << std::setprecision(10);
	std::cout << "enuc: " << enuc << "\n";
	std::cout << "1e energy: " << e1.real() << "\n";
	std::cout << "spin-Zeeman energy x: " << eSpinZeemanX.real() << "\n";
	std::cout << "spin-Zeeman energy y: " << eSpinZeemanY.real() << "\n";
	std::cout << "spin-Zeeman energy z: " << eSpinZeemanZ.real() << "\n";
	std::cout << "total 1e: " << (e1 + eSpinZeemanX + eSpinZeemanY + eSpinZeemanZ).real() << "\n";



	// ==================================== two electron matrix ===================================
	//
	// calculate Fock matrix in spinor basis
	// F = h + G + ZF
	// calc G first
	Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	Eigen::MatrixXcd K = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	Eigen::MatrixXcd Cnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	Eigen::MatrixXcd Knx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	//Eigen::MatrixXcd F = hmatBig;
	//F << hmat, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
	//	Eigen::MatrixXcd::Zero(matrixSize, matrixSize), hmat;
	for (int k=0; k<spinorSize; k++) {
		for (int l=0; l<spinorSize; l++) {
			for (int m=0; m<spinorSize; m++) {
				for (int n=0; n<spinorSize; n++) {
					C(k, l) += denMat(n, m) * fciD[k][l][m][n];
					K(k, l) -= denMat(n, m) * fciD[k][n][m][l];
					Cnx(k, l) += denMat(n, m) * fcinxD[k][l][m][n];
					Knx(k, l) -= denMat(n, m) * fcinxD[k][n][m][l];
				}
			}
		}
	}

	//std::cout << "density:\n" << denMat << "\n\n";
	//std::cout << "coulomb:\n" << C << "\n\n";
	//std::cout << "exchange:\n" << K << "\n\n";
	fockSAO += C.transpose();
	fockSAO += K.transpose();
	fnx += Cnx.transpose();
	fnx += Knx.transpose();
	const auto fockNeu = spinor.adjoint() * fockSAO * spinor;
	const auto Cspinor = spinor.adjoint() * C.transpose() * spinor;// + zeemanx + zeemany + zeemanz;
	const auto Kspinor = spinor.adjoint() * K.transpose() * spinor;// + zeemanx + zeemany + zeemanz;
	const auto fock = hmatBig + Cspinor + Kspinor + zeemanx + zeemany + zeemanz;
	//std::cout << "fock:\n" << fock << "\n\n";
	for (int i=0; i<nocc; i++) {
		std::cout << "orbital energy " << i << ": " << fock(i, i) << "\n";
		std::cout << "orbital energy " << i << ": " << fockNeu(i, i) << "\n";
	}
	std::cout << "\n";


	// calculate energy weighted density matrix W
	std::cout << "calculating energy weighted density matrix...\n";
	Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	for (int a=0; a<spinorSize; a++) {
		for (int b=0; b<spinorSize; b++) {
			//for (int c=0; c<spinorSize; c++) {
			//	for (int d=0; d<spinorSize; d++) {
			//		W(a, b) += denMat(c, a) * fockSAO(c, d) * denMat(b, d);
			//	}
			//}
			for (int i=0; i<nocc; i++) {
				//W(a, b) += epsilon[i] * std::conj(spinor(a, i)) * spinor(b, i);
				W(a, b) += fock(i, i) * std::conj(spinor(a, i)) * spinor(b, i);
			}
		}
	}

	//std::cout << "fock spinor basis:\n" << fock << "\n\n";
	//std::cout << "W:\n" << W << "\n\n";





	/*****************************************************************************
	 *                              gradient stuff                               *
	 ****************************************************************************/
	std::cout << "calculating gradient...\n";
	Eigen::MatrixXcd snxBig(spinorSize, spinorSize);
	snxBig << snx, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
		Eigen::MatrixXcd::Zero(matrixSize, matrixSize), snx;
	

	std::complex<double> gradW = 0.0;
	std::complex<double> gradF = 0.0;
	std::complex<double> gradWa = 0.0;
	std::complex<double> gradHa = 0.0;
	std::complex<double> gradCa = 0.0;
	std::complex<double> gradKa = 0.0;
	std::complex<double> gradZa = 0.0;
	std::complex<double> etot = 0.0;
	std::complex<double> ec = 0.0;
	std::complex<double> ej = 0.0;
	for (int a=0; a<spinorSize; a++) {
		for (int b=0; b<spinorSize; b++) {
			//grad += hnxBig(a, b) * denMat(b, a);
			//grad -= 2.0 * snxBig(a, b) * W(b, a);
			
			//grad += W(a, b) * snxBig(a, b);
			gradW -= W(b, a) * snxBig(a, b);
			gradF += 0.5 * denMat(a, b) * (fnx(a, b) + hnxBig(b, a) + zetotnx(b, a));
			//gradF += 0.5 * denMat(b, a) * hnxBig(a, b);

			gradHa += denMat(b, a) * hnxBig(a, b);
			gradWa -= W(b, a) * snxBig(a, b);
			gradZa += denMat(b, a) * zetotnx(a, b);

			etot += denMat(a, b) * hmatkek(a, b);
		}
	}
	///*

	std::complex<double> grad2e = 0.0;
	for (int a=0; a<spinorSize; a++) {
		for (int b=0; b<spinorSize; b++) {
			for (int c=0; c<spinorSize; c++) {
				for (int d=0; d<spinorSize; d++) {
					//grad += fcinxD[a][b][c][d] * denMat(a, b) * denMat(c, d);
					gradCa += 0.5 * denMat(b, a) * denMat(d, c) * fcinxD[a][b][c][d];
					gradKa -= 0.5 * denMat(b, a) * denMat(d, c) * fcinxD[a][d][c][b];
					
					etot += 0.5 * denMat(b, a) * denMat(d, c) * fciD[a][b][c][d];
					ec += 0.5 * denMat(b, a) * denMat(d, c) * fciD[a][b][c][d];
					etot -= 0.5 * denMat(b, a) * denMat(d, c) * fciD[a][d][c][b];
					ej -= 0.5 * denMat(b, a) * denMat(d, c) * fciD[a][d][c][b];
				}
			}
		}
	}
	etot += enuc + eSpinZeemanX + eSpinZeemanY + eSpinZeemanZ;
	//*/

	std::cout << "\ngradient (without Vnn): " << gradW + gradF << "\n";
	std::cout << "gradW:                  " << gradW << "\n";
	std::cout << "gradF+H:                " << gradF << "\n";
	std::cout << "total:                  " << gradW + gradF + enucnx[nuc][cart] << "\n";
	std::cout << "\n\nalternativer weg:\n";
	std::cout << "spin factor: " << spinFactor << "\n";
	std::cout << "gradHa:  " << gradHa << "\n";
	std::cout << "gradWa:  " << gradWa << "\n";
	std::cout << "gradCa:  " << gradCa << "\n";
	std::cout << "gradKa:  " << gradKa << "\n";
	std::cout << "gradZa:  " << gradZa << "\n";
	std::cout << "gradnuc: " << enucnx[nuc][cart] << "\n";
	std::cout << "total:   " << gradHa + gradWa + gradCa + gradKa + gradZa + enucnx[nuc][cart] << "\n";
	std::cout << "\n\ntotal SCF-energy:   " << etot << "\n";

	std::cout << "\n\nnuclear gradient: " << enucnx[nuc][cart] << "\n";
	std::cout << "vor coulomb:      " << gradHa + gradWa + gradZa + enucnx[nuc][cart] << "\n";
	std::cout << "vor exchange:     " << gradHa + gradWa + gradCa + gradZa + enucnx[nuc][cart] << "\n";
	std::cout << "total:            " << gradHa + gradWa + gradCa + gradKa + gradZa + enucnx[nuc][cart] << "\n";

	//std::cout << "\n\nW(0, 0) = " << W(0, 0) << "\n";
	//std::cout << "\n\nW:\n" << W << "\n";

	//const auto wread = readHerm("wcao");
	//Eigen::MatrixXcd wreadBig(spinorSize, spinorSize);
	//wreadBig << wread, Eigen::MatrixXcd::Zero(matrixSize, matrixSize),
	//	Eigen::MatrixXcd::Zero(matrixSize, matrixSize), wread;
	//auto wkek = spinor.adjoint() * wreadBig.transpose() * spinor;
	//std::cout << "\n\nWread:\n" << wkek << "\n";

	//std::cout << "2e gradient: " << grad2e << "\n";
	//std::cout << "total: " << grad + grad2e << "\n";
	//std::cout << "ec: "  << ec << "\n";
	//std::cout << "ej: "  << ej << "\n";
	
	// debug energy
	std::complex<double> edebug = 0.0;
	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			edebug += 0.5 * denMat(i, j) * (fockSAO(i, j) + hmatkek(i, j));
			//std::cout << 0.5 * denMat(i, j) * (fockSAO(i, j) + hmatkek(i, j)) << "\n";
		}
	}
	//std::cout << "\n\ndebug energy: " << edebug+enuc << "\n\n";

	//std::cout << std::setprecision(5);
	//std::cout << "\n\nW:\n" << W << "\n";
	//std::cout << "\n\nWaa:\n" << W(Eigen::seqN(0, matrixSize), Eigen::seqN(0, matrixSize)) << "\n";
	//std::cout << "\n\nWbb:\n" << W(Eigen::seqN(matrixSize, matrixSize), Eigen::seqN(matrixSize, matrixSize)) << "\n";
	//std::cout << "\n\nWsym:\n" << W(Eigen::seqN(0, matrixSize), Eigen::seqN(0, matrixSize)) + W(Eigen::seqN(matrixSize, matrixSize), Eigen::seqN(matrixSize, matrixSize)) << "\n";

}
