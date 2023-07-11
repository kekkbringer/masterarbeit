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

#define DEBUG 1
#define IFDBG if constexpr (DEBUG)

//#define DBOC

#define MAJOR_VERSION 0
#define MINOR_VERSION 0
#define PATCH_VERSION 1

using namespace std::complex_literals;

int main(int argc, char* argv[]) {
	const std::string cartDict[] = {"x", "y", "z"};
	constexpr bool usenum = false;

	int atomNum;
	int nocc;
	int nvirt;
	double Bx;
	double By;
	double Bz;
	double Bnorm;
	info(atomNum, nocc, nvirt, Bx, By, Bz, Bnorm);

	std::cout << ":: Basic infos\n";
	std::cout << "      number of atoms:                " << atomNum << "\n";
	std::cout << "      number of occupied spinors:     " << nocc << "\n";
	std::cout << "      number of virtual spinors:      " << nvirt << "\n";
	std::cout << "      strength of magnetic field:     " << Bnorm << "\n";
	std::cout << "      B field vector (unscaled):   " << Bx << "   " << By << "   " << Bz << "\n";

	std::vector<double> epsilon;
	const auto spinor = readSpinor(epsilon);

	const int spinorSize = nocc + nvirt;


	// CAO spinor test section
	std::cout << "starting CAO test section" << std::endl;
	makeCAOtrans();
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
	std::cout << "transmat done" << std::endl;
	Eigen::MatrixXcd transMatBIG(2*transMat.rows(), 2*transMat.cols());
	transMatBIG << transMat, Eigen::MatrixXcd::Zero(transMat.rows(), transMat.cols()),
			Eigen::MatrixXcd::Zero(transMat.rows(), transMat.cols()), transMat;
	
	const auto spinorCAO = transMatBIG * spinor;
	std::cout << "spinorCAO done: " << spinorCAO.rows() << " x " << spinorCAO.cols() << std::endl;

	const auto smatC = readHerm("smatcao");
	// reorder smatCAO
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
	const auto smatCAO = swapMat.transpose() * smatC * swapMat;
	
	Eigen::MatrixXcd smatCAOBIG(2*smatCAO.rows(), 2*smatCAO.cols());
	smatCAOBIG << smatCAO, Eigen::MatrixXcd::Zero(smatCAO.rows(), smatCAO.cols()),
			Eigen::MatrixXcd::Zero(smatCAO.rows(), smatCAO.cols()), smatCAO;
	std::cout << "smatCAOBIG: " << smatCAOBIG.rows() << " x " << smatCAOBIG.cols() << std::endl;
	const auto hofident = spinorCAO.adjoint() * smatCAOBIG * spinorCAO;
	//std::cout << "hofident:\n" << hofident << "\n\n";
	bool iden = true;
	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			if (i==j) {
				if (abs(abs(hofident(i, j))-1.0) >= 1e-12) {
					std::cout << "NOT IDENTITY!\n";
					iden = false;
					break;
				}
			} else {
				if (abs(hofident(i, j)) >= 1e-12) {
					std::cout << "NOT IDENTITY!\n";
					iden = false;
					break;
				}
			}
		}
	}
	if (iden) std::cout << "caosao trafo worked!" << std::endl;
	// end of CAO spinor test section



	const auto begin1 = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXcd A, B;
	const auto fciMO = eritrans(spinorCAO, nocc, nvirt, A, B);
	// add spinor energies to diagonal of A
	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			A( i*nvirt+a-nocc, i*nvirt+a-nocc ) += epsilon[a] - epsilon[i];
		}
	}
	

	const auto end1 = std::chrono::high_resolution_clock::now();


	const auto begin2 = std::chrono::high_resolution_clock::now();
	std::cout << "\n\ncalculating orbital rotation matrix\n\n";
	std::vector<Eigen::VectorXcd> ball;
	split1efiles(atomNum);
	for (int nuc=0; nuc<atomNum; nuc++) {
		for (int cart=0; cart<3; cart++) {
			std::cout << " :: calculating rhs...   ";
			const auto b0ai = berryRHS(nuc, cart, fciMO);
			ball.push_back(b0ai);
			std::cout << "done!\n";
		}
	}

	const auto interRHScphf = std::chrono::high_resolution_clock::now();

	std::cout << " :: solving CPHF equation...   " << std::flush;
	const auto u = cphf(A, B, ball);
	std::cout << "done!\n";
	std::cout << "      saving to disk...   " << std::flush;
	int ind=0;
	for (int nuc=0; nuc<atomNum; nuc++) {
		for (int cart=0; cart<3; cart++) {
			saveVector(u[ind], "u" + std::to_string(nuc) + "_" + std::to_string(cart));
			ind++;
		}
	}
	std::cout << "done!\n";

	std::cout << "\n\ndone calculating orbital rotation matrix\n\n";
	std::cout << std::flush;
	const auto end2 = std::chrono::high_resolution_clock::now();


	const auto begin3 = std::chrono::high_resolution_clock::now();
	// split bra ket files
	std::cout << "\nsplitting sbraket files...\n" << std::flush;
	splitBraKet(atomNum);
	std::cout << "done!\n" << std::flush;

	///* calculate actual berry curvature
	Eigen::MatrixXcd berry  = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry2 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry3 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry4 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry5 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry6 = Eigen::MatrixXcd::Zero(3*atomNum, 3*atomNum);

	Eigen::MatrixXd berryopt  = Eigen::MatrixXd::Zero(3*atomNum, 3*atomNum);

	Eigen::VectorXcd dboc(3*atomNum);

#pragma omp parallel for
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			const auto uIA = readVector("u" + std::to_string(I) + "_" + std::to_string(alpha));
			
			// bra abgeleitete matrix
			const auto snxbraIA = readMatrixTransform("b" + std::to_string(I) + cartDict[alpha]);
			//const auto snxketIA = readMatrixTransform("k" + std::to_string(I) + cartDict[alpha]);
			
			Eigen::MatrixXcd tmp1(spinorSize, spinorSize);
			tmp1 << snxbraIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
				Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbraIA;
			//const auto snxbraIAMO = spinor.adjoint() * tmp1 * spinor;
			const auto snxbraIAMOip = spinor.leftCols(nocc).adjoint() * tmp1 * spinor;

			////// test stuff und DBOC
			//Eigen::MatrixXcd tmp4(spinorSize, spinorSize);
			//tmp4 << snxketIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
			//	Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxketIA;
			//const auto snxketIAMO = spinor.adjoint() * tmp4 * spinor;


			// DBOC section
//#ifdef DBOC
//			// doppelt abgeleitete overlap matrix
//			auto braketIAIA = readMatrixTransform("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(I) + cartDict[alpha]);
//			Eigen::MatrixXcd tmpDBOC1(spinorSize, spinorSize);
//			tmpDBOC1 << braketIAIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
//				Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braketIAIA;
//			const auto braketDBOC = spinor.adjoint() * tmpDBOC1 * spinor;
//
//			dboc(3*I+alpha) = 0.0;
//			for (int i=0; i<nocc; i++) {
//				dboc(3*I+alpha) += braketDBOC(i, i);
//				for (int r=0; r<spinorSize; r++) {
//					dboc(3*I+alpha) -= abs(snxbraIAMO(i, r)) * abs(snxbraIAMO(i, r));
//				}
//				for (int a=nocc; a<spinorSize; a++) {
//					dboc(3*I+alpha) += abs(snxketIAMO(a, i) + uIA(i*nvirt+a-nocc)) * abs(snxketIAMO(a, i) + uIA(i*nvirt+a-nocc));
//				}
//			}
//#endif
			// end of DBOC section


			for (int J=0; J<=I; J++) {
				int betamax = 3;
				if (I==J) betamax = alpha;
				for (int beta=0; beta<betamax; beta++) {
					auto uJB = readVector("u" + std::to_string(J) + "_" + std::to_string(beta));
					
					// bra abgeleitete matrix
					//auto snxbraJB = readMatrixTransform("b" + std::to_string(J) + cartDict[beta]);
					auto snxketJB = readMatrixTransform("k" + std::to_string(J) + cartDict[beta]);
					
					Eigen::MatrixXcd tmp2(spinorSize, spinorSize);
					tmp2 << snxketJB, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
						Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxketJB;
					//const auto snxketJBMO = spinor.adjoint() * tmp2 * spinor;
					const auto snxketJBMOip = spinor.adjoint() * tmp2 * spinor.leftCols(nocc);

					//Eigen::MatrixXcd tmpa(spinorSize, spinorSize);
					//tmpa << snxbraJB, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
					//	Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbraJB;
					//const auto snxbraJBMO = spinor.adjoint() * tmpa * spinor;
					

					// doppelt abgeleitete overlap matrix
					auto braketA = readMatrixTransform("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta]);

					Eigen::MatrixXcd tmp3(spinorSize, spinorSize);
					tmp3 << braketA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
						Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braketA;
					//const auto braketMO = spinor.adjoint() * tmp3 * spinor;


					for (int i=0; i<nocc; i++) {
						//berry(3*I+alpha, 3*J+beta) += braketMO(i, i);
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
					
					berry(3*I+alpha, 3*J+beta) *= -2;
					berry2(3*I+alpha, 3*J+beta) *= -2;
					berry3(3*I+alpha, 3*J+beta) *= -2;
					berry4(3*I+alpha, 3*J+beta) *= -2;
					berry5(3*I+alpha, 3*J+beta) *= -2;
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
	const auto end3 = std::chrono::high_resolution_clock::now();

	std::cout << std::fixed;
	std::cout << std::setprecision(10);
	//std::cout << "\n\n\nBerry-curvature term 1:\n" <<  berry.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 2:\n" << berry2.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 3:\n" << berry3.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 4:\n" << berry4.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 5:\n" << berry5.imag() << "\n\n";
	std::cout << "\n\n\nBerry-curvature total:\n" <<  berry6.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature sym:\n" << 0.5*(berry + berry.transpose()).imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature antisym:\n" << 0.5*(berry - berry.transpose()).imag() << "\n\n";
	
	std::cout << "================================================================================\n";

	// calculate shielding charges
	Eigen::MatrixXd chargeFluct = Eigen::MatrixXd::Zero(atomNum, atomNum);
	for (int I=0; I<atomNum; I++) {
		for (int J=0; J<atomNum; J++) {
			// get IJ block from total berry curvature
			const Eigen::MatrixXd oij = berry6.block<3,3>(3*I, 3*J).imag();
			//std::cout << "IJ = " << I << "  " << J << "\n";
			//std::cout << oij << "\n";

			// get antisymmetric part
			const auto oijas = 0.5*(oij-oij.transpose());
			//std::cout << "\n" << oijas << "\n";

			//const Eigen::Vector3d omegaij(oijas(2, 1), oijas(0, 2), oijas(1, 0));
			//std::cout << "\n" << omegaij << "\n";

			chargeFluct(I, J) = oijas(2, 1) * Bx
					  + oijas(0, 2) * By
					  + oijas(1, 0) * Bz;
			//chargeFluct(I, J) /= Bnorm * Bnorm;
			chargeFluct(I, J) /= Bx*Bx + By*By + Bz*Bz;
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
		//std::cout << " " << I+1 << atoms[I] << "\t" << q << "\t" << chargeOf(atoms[I]) << "\t\t" << q+chargeOf(atoms[I]) << "\n";
		printf(" %u%s\t%10.7f\t%10.7f\t%10.7f\n", I+1, atoms[I].c_str(), q, (double)chargeOf(atoms[I]), q+chargeOf(atoms[I]));
		esum += q;
		nsum += chargeOf(atoms[I]);
	}
	std::cout << "-----------------------------------------------------\n";
	//std::cout << " sum\t" << esum << "\t" << nsum << "\t" << esum+nsum << "\n";
	printf(" sum\t%10.7f\t%10.7f\t%10.7f\n", esum, nsum, esum+nsum);
	//*/
	
#ifdef DBOC
	// print DBOC
	double dboctot = 0.0;
	//std::cout << "\n\nDBOC:\n";
	for (int I=0; I<atomNum; I++) {
		//std::cout << "  atom " << I+1 << ":\n";
		//std::cout << "    x: " << dboc(3*I + 0)/(2*massOf(atoms[I])*1822.8885291649) << "H   = " << 219474.6 * dboc(3*I + 0)/(2*massOf(atoms[I])*1822.8885291649) << "cm-1\n";
		//std::cout << "    y: " << dboc(3*I + 1)/(2*massOf(atoms[I])*1822.8885291649) << "H   = " << 219474.6 * dboc(3*I + 1)/(2*massOf(atoms[I])*1822.8885291649) << "cm-1\n";
		//std::cout << "    z: " << dboc(3*I + 2)/(2*massOf(atoms[I])*1822.8885291649) << "H   = " << 219474.6 * dboc(3*I + 2)/(2*massOf(atoms[I])*1822.8885291649) << "cm-1\n";
		dboctot += 	  dboc(3*I + 0).real()/(2*massOf(atoms[I])*1822.8885291649)
				+ dboc(3*I + 1).real()/(2*massOf(atoms[I])*1822.8885291649)
				+ dboc(3*I + 2).real()/(2*massOf(atoms[I])*1822.8885291649);
	}
	std::cout << "\n\ntotal DBOC = " << dboctot << " H  = " << dboctot*219474.6 << " cm-1";
	std::cout << "\n\n";
#endif
	
	
	// time stats
	const auto elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1-begin1);
	const auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(interRHScphf-begin2);
	const auto elapsedi = std::chrono::duration_cast<std::chrono::milliseconds>(end2-interRHScphf);
	const auto elapsed3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3-begin3);
	std::cout << "\n\n =================== time stats ===================\n";
	printf("   Calculate electronic Hessian: %.3fs\n", elapsed1.count()*1e-3);
	printf("   Calculate RHS:                %.3fs\n", elapsed2.count()*1e-3);
	printf("   Solve CPHF:                   %.3fs\n", elapsedi.count()*1e-3);
	printf("   Calculate Berry-Curvature:    %.3fs\n", elapsed3.count()*1e-3);
	printf("   --------------------------------------------------------\n");
	printf("   Total:                        %.3fs\n\n", (elapsed1.count()+elapsed2.count()+elapsed3.count()) * 1e-3);
	
	
	std::cout << std::endl;
	return 0;
}
