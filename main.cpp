#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <complex>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <algorithm>

#include "misc.hpp"
#include "cphf.hpp"
#include "read_spinor.hpp"
#include "read_herm.hpp"
#include "read_fourcenter.hpp"
#include "calc_stabmat.hpp"
#include "fci_grad.hpp"
#include "berry_rhs.hpp"
#include "misc.hpp"
#include "eritrans.hpp"
#include "matrixio.hpp"
#include "restart.hpp"

//#define DBOC

using namespace std::complex_literals;

/* 
 * The main function that drives the overall computation of the Berry curvature
 * and Berry charges. Most calculations are done in their own function, some are
 * still left here. Still the structure of the main function can largely be
 * divided into four main steps:
 * 	1. Calculation of the electronic hessian
 * 	2. Calculation of the CPHF right-hand side vectors
 * 	3. Solving the CPHF equations
 * 	4. Calculating the Berry curvature and related quantities
 *
 * Each section will be preluded by a small explanation of its most important
 * features and functionalities.
 *
 */
int main(int argc, char* argv[]) {
	setRestart(0);

	/**********************************************************************
	 *                         Electronic Hessian                         *
	 **********************************************************************
	 * In this section basic information about the molecular geometry and
	 * basis set are read. In addition the CAO basis spinor coefficients
	 * are calculated, which are needed for the contraction of fourcenter
	 * quantities.
	 * Next, the electronic Hessian
	 *
	 *	    /  A   B  \
	 * 	H = |         |
	 * 	    \  A*  B* /
	 *
	 * is computed by a batched transformation of fourcenter LAOs. During
	 * the transformation integrals of the type (ai|lk) are saved on the
	 * heap as they are later reused to compute the rhs vectors for CPHF.
	 *
	 * The fourcenter LAOs needed here are read from files 'fourcenter.r'
	 * and 'fourcenter.i'.
	 */
	const auto begin1 = std::chrono::high_resolution_clock::now();

	const std::string cartDict[] = {"x", "y", "z"};

	// get basis info about the molecule, basis and magentic field
	int atomNum, nocc, nvirt;
	double Bx, By, Bz, Bnorm;
	info(atomNum, nocc, nvirt, Bx, By, Bz, Bnorm);
	std::cout << ":: Basic infos\n";
	std::cout << "      number of atoms:                " << atomNum << "\n";
	std::cout << "      number of occupied spinors:     " << nocc << "\n";
	std::cout << "      number of virtual spinors:      " << nvirt << "\n";
	std::cout << "      strength of magnetic field:     " << Bnorm << "\n";
	std::cout << "      B field vector (unscaled):   " << Bx << "   " << By << "   " << Bz << "\n";

	// read spinor coefficients and spinor energies
	std::vector<double> epsilon;
	const auto spinor = readSpinor(epsilon);
	const int spinorSize = nocc + nvirt;

	// Compute CAO spinor coefficients
	makeCAOtrans(); // calculate reordering of AOs and transformation matrix for CAO <-> SAO
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
	//std::cout << "transmat done" << std::endl;
	Eigen::MatrixXcd transMatBIG(2*transMat.rows(), 2*transMat.cols());
	transMatBIG << transMat, Eigen::MatrixXcd::Zero(transMat.rows(), transMat.cols()),
			Eigen::MatrixXcd::Zero(transMat.rows(), transMat.cols()), transMat;
	
	const auto spinorCAO = transMatBIG * spinor;
	std::cout << "CAO spinor done: " << spinorCAO.rows() << " x " << spinorCAO.cols() << "\n\n" << std::endl;
	const int ncao = spinorCAO.rows()/2;

	// compute electronic hessian
	Eigen::MatrixXcd A, B;
	const auto fciMO = eritrans(spinorCAO, nocc, nvirt, A, B);
	// add spinor energies to diagonal of A
	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			A( i*nvirt+a-nocc, i*nvirt+a-nocc ) += epsilon[a] - epsilon[i];
		}
	}

	const auto end1 = std::chrono::high_resolution_clock::now();
	setRestart(1);
	// end of electronic Hessian section


	/**********************************************************************
	 *                              CPHF rhs                              *
	 **********************************************************************
	 * This section is concerned with the calculation of the right-hand
	 * side vectors of the Coupled Perturbed Hartree-Fock equations.
	 *
	 * 	b_ai = F'_ai + S'*eps + S'(ai||lk)
	 *
	 * The function 'berryRHS' drives their calculation, 1e contributions
	 * are directly computed in this function. 2e contributions to the der-
	 * ivative Fock matrix are computed in TURBOMOLE and read from files
	 * 'gjfock.r' and 'gjfock.i' (Coulomb derivatives) and 'gxfock' (Ex-
	 * change derivatives). 2e contributions contracted with the derivative
	 * overlap matrix are calculated with the spinor basis fourcenter inte-
	 * grals computed in the previous section. The rhs vectors are written
	 * to disk on files 'b0ai...' and simultaniously collected and kept in
	 * memory.
	 *
	 * Derivative overlap matrices are read from files 'sbra.r/.i' and
	 * 'sket.r/.i', where only the bra- or ket- basis function is a deri-
	 * vative. The derivative core Hamilton is read from 'hgrad.r/.i',
	 * which already includes all 1e contributions, including the Hellman-
	 * Feynman contribtution. Only wavefunction depended 1e contributions
	 * like the spin-Zeeman are not yet included.
	 */

	const auto begin2 = std::chrono::high_resolution_clock::now();
	std::vector<Eigen::VectorXcd> ball;
	splitExchange(atomNum, ncao);
	split1efiles(atomNum, ncao);
	std::cout << "\n\n :: calculating CPHF rhs vectors\n\n" << std::flush;
	for (int nuc=0; nuc<atomNum; nuc++) {
		for (int cart=0; cart<3; cart++) {
			const auto beginrhs = std::chrono::high_resolution_clock::now();
			const auto b0ai = berryRHS(nuc, cart, fciMO);
			ball.push_back(b0ai);
			const auto endrhs = std::chrono::high_resolution_clock::now();
			const auto elapsedrhs = std::chrono::duration_cast<std::chrono::milliseconds>(endrhs-beginrhs);
			std::cout << "after " << elapsedrhs.count()*1e-3 << " s\n\n" << std::endl;
		}
	}
	const auto interRHScphf = std::chrono::high_resolution_clock::now();
	setRestart(2);
	// end of rhs section

	
	/**********************************************************************
	 *                             CPHF solver                            *
	 **********************************************************************
	 * In this section the CPHF equations are solved. The method of choice
	 * is the Davidson algorithm, that solveds the CPHF equations for all
	 * rhs vectors simultaneously by projecting the linear system of equa-
	 * tions to a subspace of increasing dimension until the subspace spans
	 * the same space spanned by the full SLE up to a certain threshold.
	 *
	 * 	/  A  B  \ / b  \     / U  \
	 * 	|        | |    |  =  |    |
	 * 	\  A* B* / \ b* /     \ U* /
	 *
	 * The solution vectors U are written to disk in files 'u...' and read
	 * again, later.
	 */
	std::cout << " :: solving CPHF equation...   " << std::flush;
	const auto u = cphf(A, B, ball); // solving CPHF
	std::cout << "\tsaving to disk...   " << std::flush;
	int ind=0;
	for (int nuc=0; nuc<atomNum; nuc++) {
		for (int cart=0; cart<3; cart++) {
			saveVector(u[ind], "u" + std::to_string(nuc) + "_" + std::to_string(cart));
			ind++;
		}
	}
	std::cout << "done!\n";
	std::cout << std::flush;
	const auto end2 = std::chrono::high_resolution_clock::now();
	setRestart(3);
	// end of CPHF section


	/**********************************************************************
	 *                    Berry curvature and charges                     *
	 **********************************************************************
	 * This section drives the calculation of the Berry curavture and
	 * charges. 
	 * In the first step the overlap matrices <mu'|nu>, <mu|nu'>
	 * and <mu’|nu’> are transformed to spinor basis and stored on disk.
	 * For the second step the lower triangle of the Berry tensor is
	 * calculated acoording to eq. 2.86 (master's thesis Steinmetz).
	 * Here, all required quantities can be obtained from disk: The deri-
	 * vative overlap matrices in spinor basis and the CPHF solution vec-
	 * tors U are all read as needed. The indiviual terms in eq. 2.86
	 * are all treated differently, their real and imaginary parts are
	 * calculated even though only the imaginary part is needed for the
	 * Berry tensor itself.
	 */
	const auto begin3 = std::chrono::high_resolution_clock::now();
	// split bra ket files
	std::cout << "\nsplitting sbraket files..." << std::flush;
	splitBraKet(atomNum);
	std::cout << "   done!\n" << std::flush;

	// transform derivative overlap matrices to spinor basis
	const auto begintrafo = std::chrono::high_resolution_clock::now();
	std::cout << "transforming integrals to spinor basis..." << std::flush;
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			// bra derivative
			const auto snxbraIA = readMatrixTransform("b" + std::to_string(I) + cartDict[alpha]);
			Eigen::MatrixXcd tmp1(spinorSize, spinorSize);
			tmp1 << snxbraIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
				Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbraIA;
			const auto snxbraIAMOip = spinor.leftCols(nocc).adjoint() * tmp1 * spinor;
			std::string filename = "snxbraip" + std::to_string(I) + "_" + std::to_string(alpha);
			saveMatrix(snxbraIAMOip, filename);

			// ket derivative
			const auto snxketIA = readMatrixTransform("k" + std::to_string(I) + cartDict[alpha]);
			Eigen::MatrixXcd tmp2(spinorSize, spinorSize);
			tmp2 << snxketIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
				Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxketIA;
			const auto snxketIAMOip = spinor.adjoint() * tmp2 * spinor.leftCols(nocc);
			filename = "snxketip" + std::to_string(I) + "_" + std::to_string(alpha);
			saveMatrix(snxketIAMOip, filename);
		}
	}
	const auto endtrafo = std::chrono::high_resolution_clock::now();
	const auto elapsedtrafo = std::chrono::duration_cast<std::chrono::milliseconds>(endtrafo-begintrafo);
	std::cout << "   done after " << elapsedtrafo.count()*1e-3 << " s" << std::endl;



	// calculate actual berry curvature
	Eigen::MatrixXcd berry  = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry2 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry3 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry4 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry5 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry6 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);

#pragma omp parallel for
	std::cout << "\ncalculating Berry curvature tensor..." << std::endl;
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			// U vector I alpha
			const auto uIA = readVector("u" + std::to_string(I) + "_" + std::to_string(alpha));
			
			// bra derivative overlap matrix
			const std::string filename = "snxbraip" + std::to_string(I) + "_" + std::to_string(alpha);
			const auto snxbraIAMOip = loadMatrix(filename);

			for (int J=0; J<=I; J++) {
				int betamax = 3;
				if (I==J) betamax = alpha;
				for (int beta=0; beta<betamax; beta++) {
					// U vector J beta
					auto uJB = readVector("u" + std::to_string(J) + "_" + std::to_string(beta));
					
					// ket derivative overlap matrix
					const std::string filename2 = "snxketip" + std::to_string(J) + "_" + std::to_string(beta);
					const auto snxketJBMOip = loadMatrix(filename2);

					// bra and ket derivative overlap matrix
					auto braketA = readMatrixTransform("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta]);

					// spin up/up and down/down combinations (inefficient)
					Eigen::MatrixXcd tmp3(spinorSize, spinorSize);
					tmp3 << braketA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
						Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braketA;

					for (int i=0; i<nocc; i++) {
						berry(3*I+alpha, 3*J+beta) += (spinor.col(i).adjoint() * tmp3 * spinor.col(i))(0, 0);
						
						for (int a=nocc; a<spinorSize; a++) {
							berry2(3*I+alpha, 3*J+beta) += snxketJBMOip(a, i) * std::conj(uIA(i*nvirt+a-nocc));
							berry3(3*I+alpha, 3*J+beta) += snxbraIAMOip(i, a) * uJB(i*nvirt+a-nocc);
							berry4(3*I+alpha, 3*J+beta) += uIA(i*nvirt+a-nocc + nocc*nvirt) * uJB(i*nvirt+a-nocc);
						}
						for (int j=0; j<nocc; j++) {
							berry5(3*I+alpha, 3*J+beta) -= snxbraIAMOip(i, j) * snxketJBMOip(j, i);
						}
					}
					berry6(3*I+alpha, 3*J+beta) = berry(3*I+alpha, 3*J+beta)
									+ berry2(3*I+alpha, 3*J+beta)
									+ berry3(3*I+alpha, 3*J+beta)
									+ berry4(3*I+alpha, 3*J+beta)
									+ berry5(3*I+alpha, 3*J+beta);
					berry6(3*J+beta, 3*I+alpha) = std::conj(berry6(3*I+alpha, 3*J+beta));
				}
			}
		}
	}
	berry *= -2;	berry2 *= -2;	berry3 *= -2;
	berry4 *= -2;	berry5 *= -2;	berry6 *= -2;
	const auto end3 = std::chrono::high_resolution_clock::now();

	std::cout << std::fixed;
	std::cout << std::setprecision(10);
	//std::cout << "\n\n\nBerry-curvature term 1:\n" <<  berry.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 2:\n" << berry2.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 3:\n" << berry3.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 4:\n" << berry4.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 5:\n" << berry5.imag() << "\n\n";
	std::cout << "\n\n\nBerry-curvature total:\n" <<  berry6.imag() << "\n\n";
	std::cout << "================================================================================\n";
	saveBerry(berry6);

	// calculate shielding charges
	Eigen::MatrixXd chargeFluct = Eigen::MatrixXd::Zero(atomNum, atomNum);
	for (int I=0; I<atomNum; I++) {
		for (int J=0; J<atomNum; J++) {
			// get IJ block from total berry curvature
			const Eigen::MatrixXd oij = berry6.block<3,3>(3*I, 3*J).imag();

			// get antisymmetric part
			const auto oijas = 0.5*(oij-oij.transpose());

			chargeFluct(I, J) = oijas(2, 1) * Bx
					  + oijas(0, 2) * By
					  + oijas(1, 0) * Bz;
			chargeFluct(I, J) /= Bx*Bx + By*By + Bz*Bz; // caution is needed, this only really works for magnetic fields in z direction!
			//chargeFluct(I, J) *= 2.35051756758e5; // TESLA, ÄNDERE DAS NOCH UNBEDINGT, DU IDIOT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
	}
	std::cout << "\n\ncharge fluctuation:\n" << chargeFluct << "\n\n";


	// kurz ins coord file nach den Atomen gucken... :D
	std::ifstream coord("coord");
	//std::string line, word;
	std::string word;
	std::vector<std::string> atoms;
	getline(coord, line); //$coord
	for (int i=0; i<atomNum; i++) {
		getline(coord, line);
		std::istringstream iss(line);
		iss >> word; // x
		iss >> word; // y
		iss >> word; // z
		iss >> word; // atom
		atoms.push_back(word);
	}
	coord.close();


	// calculate resulting partial charges
	std::cout << "\npartial charges:\n";
	std::cout << "Atom\telec. charge \tnuc. charge\ttotal\n";
	std::vector<double> partialCharges;
	double esum=0.0, nsum=0.0;
	for (int I=0; I<atomNum; I++) {
		double q = 0.0;
		for (int J=0; J<atomNum; J++) {
			q += chargeFluct(I, J);
		}
		partialCharges.push_back(q);
		printf(" %u%s\t% 10.7f\t% 10.7f\t% 10.7f\n", I+1, atoms[I].c_str(), q, (double)chargeOf(atoms[I]), q+chargeOf(atoms[I]));
		esum += q;
		nsum += chargeOf(atoms[I]);
	}
	std::cout << "-----------------------------------------------------\n";
	printf(" sum\t% 10.7f\t% 10.7f\t% 10.7f\n", esum, nsum, esum+nsum);
	
	
	// time stats
	const auto elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1-begin1);
	const auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(interRHScphf-begin2);
	const auto elapsedi = std::chrono::duration_cast<std::chrono::milliseconds>(end2-interRHScphf);
	const auto elapsed3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3-begin3);
	std::cout << "\n\n =================== time stats ===================\n";
	printf("   Calculate electronic Hessian:   %.3f s\n", elapsed1.count()*1e-3);
	printf("   Calculate RHS:                  %.3f s\n", elapsed2.count()*1e-3);
	printf("   Solve CPHF:                     %.3f s\n", elapsedi.count()*1e-3);
	printf("   Calculate Berry-Curvature:      %.3f s\n", elapsed3.count()*1e-3);
	printf(" --------------------------------------------------\n");
	printf("   Total:                          %.3f s\n\n", (elapsed1.count()+elapsed2.count()+elapsed3.count()) * 1e-3);

	// delete temporary scratch files
	deleteTmpFiles(atomNum);
	std::cout << std::endl;
	// end of Berry section


	return 0;
}
