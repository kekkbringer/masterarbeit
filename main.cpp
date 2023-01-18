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
#include "berry_rhs.hpp"
#include "misc.hpp"

#define DEBUG 1
#define IFDBG if constexpr (DEBUG)

#define MAJOR_VERSION 0
#define MINOR_VERSION 0
#define PATCH_VERSION 1

using namespace std::complex_literals;

int main(int argc, char* argv[]) {
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
	std::cout << "      B field vector:   " << Bx << "   " << By << "   " << Bz << "\n";

	std::vector<double> epsilon;
	const auto spinor = readSpinor(epsilon);

	const int spinorSize = nocc + nvirt;
	
	Eigen::MatrixXcd A;
	Eigen::MatrixXcd B;
	calcStabmat(A, B);

	std::cout << "\n\n\n\n\n\ncalculating orbital rotation matrix\n\n\n\n\n\n\n";
	const std::string cartDict[] = {"x", "y", "z"};
	for (int nuc=0; nuc<atomNum; nuc++) {
		for (int cart=0; cart<3; cart++) {
			std::cout << " :: calculating rhs...   ";
			const auto b0ai = berryRHS(nuc, cart);
			std::cout << "done!\n";

			std::cout << " :: solving CPHF equation...   ";
			const auto u = cphf(A, B, b0ai);
			//std::cout << "\nU-Vector:\n" << u << "\n\n";
			std::cout << "done!\n";
			std::cout << "      saving to disk...   ";
			saveVector(u, "u" + std::to_string(nuc) + "_" + std::to_string(cart));
			std::cout << "done!\n";

			// TBD
			//Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
			//for (int i=0; i<nocc; i++) {
			//	for (int a=nocc; a<spinorSize; a++) {
			//		U(a, i) = u(i*nvirt+a-nocc);
			//	}
			//}
			//std::cout << "\nU:\n" << U << "\n\n";

			//auto snxbra = readHerm("b" + std::to_string(nuc) + cartDict[cart]);
			//auto snxket = readHerm("k" + std::to_string(nuc) + cartDict[cart]);
		}
	}

	std::cout << "\n\n\n\n\n\ndone calculating orbital rotation matrix\n\n\n\n\n\n\n";

	// split bra ket files
	splitBraKet(atomNum);

	///* calculate actual berry curvature
	Eigen::MatrixXcd berry(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry2(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry3(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry4(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry5(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry5bb(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry5kk(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry6(3*atomNum, 3*atomNum);
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			Eigen::MatrixXcd uNumIA;
			try {
				auto cplusIA  = readNumSpinor("TEST5/" + std::to_string(3*I+alpha+1) + "IPlus.out");
				auto cminusIA = readNumSpinor("TEST5/" + std::to_string(3*I+alpha+1) + "IMinus.out");
				auto cnxNumIA = (cplusIA - cminusIA) / 2e-3;
				uNumIA = spinor.inverse() * cnxNumIA;
			}
			catch (...) {}

			//std::cout << "U:\n" << std::fixed << std::setprecision(4) << uNumIA << "\n\n";

			auto uIA = readVector("u" + std::to_string(I) + "_" + std::to_string(alpha));
			//uIA *= -1.0;
			//std::cout << "u-vector read:\n" << uIA << "\n\n";

			//auto snxbraIA = readHerm("b" + std::to_string(I) + cartDict[alpha]);
			//Eigen::MatrixXcd tmp1(spinorSize, spinorSize);
			//tmp1 << snxbraIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
			//	Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbraIA;
			//const auto snxbraIAMO = spinor.adjoint() * tmp1 * spinor;
			
			// bra abgeleitete matrix
			auto snxbraIA = readHerm("b" + std::to_string(I) + cartDict[alpha]);
			auto snxketIA = readHerm("k" + std::to_string(I) + cartDict[alpha]);
			//auto snxbraIA = readMatrix("b" + std::to_string(I) + cartDict[alpha]);
			//auto snxketIA = readMatrix("k" + std::to_string(I) + cartDict[alpha]);
			// switch lower triangle of snxbra and snxket
			auto tmpbraIA = snxbraIA;
			for (int i=0; i<spinorSize/2; i++) {
				for (int j=0; j<i; j++) {
					snxbraIA(i, j) = snxketIA(i, j);
					snxketIA(i, j) = tmpbraIA(i, j);
				}
			}
			
			Eigen::MatrixXcd tmp1(spinorSize, spinorSize);
			tmp1 << snxbraIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
				Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbraIA;
			const auto snxbraIAMO = spinor.adjoint() * tmp1 * spinor;

			////std::cout << "snxBraIAMO:\n" << std::fixed << std::setprecision(4) << snxbraIAMO << "\n\n";

			//// test stuff
			Eigen::MatrixXcd tmp4(spinorSize, spinorSize);
			tmp4 << snxketIA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
				Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxketIA;
			const auto snxketIAMO = spinor.adjoint() * tmp4 * spinor;


			for (int J=0; J<atomNum; J++) {
				for (int beta=0; beta<3; beta++) {
					//std::cout << "Ia " << I << "_" << alpha << "   Jb " << J << "_" << beta << "\n";

					Eigen::MatrixXcd uNumJB;
					try {
						auto cplusJB  = readNumSpinor("TEST5/" + std::to_string(3*J+beta+1) + "IPlus.out");
						auto cminusJB = readNumSpinor("TEST5/" + std::to_string(3*J+beta+1) + "IMinus.out");
						auto cnxNumJB = (cplusJB - cminusJB) / 2e-3;
						uNumJB = spinor.inverse() * cnxNumJB;
					}
					catch (...) {}

					auto uJB = readVector("u" + std::to_string(J) + "_" + std::to_string(beta));
					//uJB *= -1.0;

					//auto snxketJB = readHerm("k" + std::to_string(J) + cartDict[beta]);
					
					// bra abgeleitete matrix
					auto snxbraJB = readHerm("b" + std::to_string(J) + cartDict[beta]);
					auto snxketJB = readHerm("k" + std::to_string(J) + cartDict[beta]);
					//auto snxbraJB = readMatrix("b" + std::to_string(J) + cartDict[beta]);
					//auto snxketJB = readMatrix("k" + std::to_string(J) + cartDict[beta]);
					// switch lower triangle of snxbra and snxket
					const auto tmpbraJB = snxbraJB;
					for (int i=0; i<spinorSize/2; i++) {
						for (int j=0; j<i; j++) {
							snxbraJB(i, j) = snxketJB(i, j);
							snxketJB(i, j) = tmpbraJB(i, j);
						}
					}
					
					Eigen::MatrixXcd tmp2(spinorSize, spinorSize);
					tmp2 << snxketJB, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
						Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxketJB;
					const auto snxketJBMO = spinor.adjoint() * tmp2 * spinor;

					Eigen::MatrixXcd tmpa(spinorSize, spinorSize);
					tmpa << snxbraJB, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
						Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbraJB;
					const auto snxbraJBMO = spinor.adjoint() * tmpa * spinor;
					

					//const auto braket = readHerm("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta]);
					//Eigen::MatrixXcd tmp3(spinorSize, spinorSize);
					//tmp3 << braket, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
					//	Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braket;
					//const auto braketMO = spinor.adjoint() * tmp3 * spinor;
					
					// doppelt abgeleitete overlap matrix
					auto braketA = readHerm("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta]);
					auto braketB = readHerm("bk" + std::to_string(J) + cartDict[beta] + std::to_string(I) + cartDict[alpha]);
					// switch lower triangle
					auto tmpbraket = braketA;
					for (int i=0; i<spinorSize/2; i++) {
						for (int j=0; j<i; j++) {
							braketA(i, j) = braketB(i, j);
							braketB(i, j) = tmpbraket(i, j);
						}
					}
					Eigen::MatrixXcd tmp3(spinorSize, spinorSize);
					tmp3 << braketA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
						Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braketA;
					const auto braketMO = spinor.adjoint() * tmp3 * spinor;


					berry(3*I+alpha, 3*J+beta) = 0.0;
					berry2(3*I+alpha, 3*J+beta) = 0.0;
					berry3(3*I+alpha, 3*J+beta) = 0.0;
					berry4(3*I+alpha, 3*J+beta) = 0.0;
					berry5(3*I+alpha, 3*J+beta) = 0.0;
					//berry5bb(3*I+alpha, 3*J+beta) = 0.0;
					//berry5kk(3*I+alpha, 3*J+beta) = 0.0;
					berry6(3*I+alpha, 3*J+beta) = 0.0;
					for (int i=0; i<nocc; i++) {
						berry(3*I+alpha, 3*J+beta) += braketMO(i, i);

						for (int a=nocc; a<spinorSize; a++) {
							if (argc<=1) berry2(3*I+alpha, 3*J+beta) += snxketJBMO(a, i) * std::conj(uNumIA(a, i));
							if (argc<=1) berry3(3*I+alpha, 3*J+beta) += snxbraIAMO(i, a) * uNumJB(a, i);
							if (argc<=1) berry4(3*I+alpha, 3*J+beta) += std::conj(uNumIA(a, i)) * uNumJB(a, i);
							
							if (argc>1) berry2(3*I+alpha, 3*J+beta) += snxketJBMO(a, i) * std::conj(uIA(i*nvirt+a-nocc));
							if (argc>1) berry3(3*I+alpha, 3*J+beta) += snxbraIAMO(i, a) * uJB(i*nvirt+a-nocc);
							if (argc>1) berry4(3*I+alpha, 3*J+beta) += uIA(i*nvirt+a-nocc + nocc*nvirt) * uJB(i*nvirt+a-nocc);
							
							// alternative für oben (nicht ganz)
							//berry(3*I+alpha, 3*J+beta) += std::conj(snxketIAMO(a, i) + uIA(i*nvirt+a-nocc)) * (snxketJBMO(a, i) + uJB(i*nvirt+a-nocc));
						}
						for (int j=0; j<nocc; j++) {
							berry5(3*I+alpha, 3*J+beta) -= snxbraIAMO(i, j) * snxketJBMO(j, i);
							//berry5bb(3*I+alpha, 3*J+beta) -= snxbraIAMO(i, j) * std::conj(snxbraJBMO(i, j));
							//berry5kk(3*I+alpha, 3*J+beta) -= std::conj(snxketIAMO(j, i)) * snxketJBMO(j, i);
							//std::cout << "term 5 -= " << snxbraIAMO(i, j) << "  *  " <<  snxketJBMO(j, i) << "  =  " << snxbraIAMO(i, j) * snxketJBMO(j, i) << "\n";
						}
						// alternative (wegen oben)
						//for (int r=0; r<spinorSize; r++) {
						//	berry(3*I+alpha, 3*J+beta) -= snxbraIAMO(i, r) * snxketJBMO(r, i);
						//}
					}

					//std::cout << "Term3 = " << berry3(3*I+alpha, 3*J+beta) << "\n\n";
					//std::cout << "\n\n";

					
					//berry(3*I+alpha, 3*J+beta) *= -2;
					//berry2(3*I+alpha, 3*J+beta) *= -2;
					//berry3(3*I+alpha, 3*J+beta) *= -2;
					//berry4(3*I+alpha, 3*J+beta) *= -2;
					//berry5(3*I+alpha, 3*J+beta) *= -2;
					berry6(3*I+alpha, 3*J+beta) = berry(3*I+alpha, 3*J+beta) + berry2(3*I+alpha, 3*J+beta) + berry3(3*I+alpha, 3*J+beta) + berry4(3*I+alpha, 3*J+beta) + berry5(3*I+alpha, 3*J+beta);
				}
			}
		}
	}
	std::cout << std::fixed;
	std::cout << std::setprecision(10);
	std::cout << "\n\n\nBerry-curvature term 1:\n" << 2.0*berry.imag() << "\n\n";
	std::cout << "\n\n\nBerry-curvature term 2:\n" << 2.0*berry2.imag() << "\n\n";
	std::cout << "\n\n\nBerry-curvature term 3:\n" << 2.0*berry3.imag() << "\n\n";
	std::cout << "\n\n\nBerry-curvature term 4:\n" << 2.0*berry4.imag() << "\n\n";
	std::cout << "\n\n\nBerry-curvature term 5:\n" << 2.0*berry5.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 5 bra bra:\n" << 2.0*berry5bb.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature term 5 ket ket:\n" << 2.0*berry5kk.imag() << "\n\n";
	std::cout << "\n\n\nBerry-curvature total:\n" <<  2.0*berry6.imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature sym:\n" << 0.5*(berry + berry.transpose()).imag() << "\n\n";
	//std::cout << "\n\n\nBerry-curvature antisym:\n" << 0.5*(berry - berry.transpose()).imag() << "\n\n";
	//*/
	


	/*
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			for (int J=0; J<atomNum; J++) {
				for (int beta=0; beta<3; beta++) {
					std::cout << I << " " << alpha << "     " << J << " " << beta << "\n";
					// doppelt abgeleitete overlap matrix
					auto braketA = readHerm("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta]);
					auto braketB = readHerm("bk" + std::to_string(J) + cartDict[beta] + std::to_string(I) + cartDict[alpha]);
					// switch lower triangle
					auto tmpbraket = braketA;
					for (int i=0; i<spinorSize/2; i++) {
						for (int j=0; j<i; j++) {
							braketA(i, j) = braketB(i, j);
							braketB(i, j) = tmpbraket(i, j);
						}
					}
					Eigen::MatrixXcd tmp3(spinorSize, spinorSize);
					tmp3 << braketA, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
						Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braketA;
					const auto braketMO = spinor.adjoint() * tmp3 * spinor;

					std::cout << "AO small real:\n" << std::fixed << std::setprecision(7) << braketA.real() << "\n\n";
					std::cout << "AO small real:\n" << std::fixed << std::setprecision(7) << braketB.real() << "\n\n";
					std::cout << "AO small imag:\n" << std::fixed << std::setprecision(7) << braketA.imag() << "\n\n";
					std::cout << "AO small imag:\n" << std::fixed << std::setprecision(7) << braketB.imag() << "\n\n";
					//std::cout << "MO real:\n" << std::fixed << std::setprecision(7) << braketMO.real() << "\n\n";
					//std::cout << "MO imag:\n" << std::fixed << std::setprecision(7) << braketMO.imag() << "\n\n";

					std::cout << "\n\n\n";
				}
			}
		}
	}
	//*/
	
	/*
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			std::cout << I << " " << alpha << "\n";


			auto snxbraIA = readMatrix("b" + std::to_string(I) + cartDict[alpha]);
			std::cout << std::fixed << std::setprecision(6) << "bra real:\n" << snxbraIA.real() << "\n";
			std::cout << std::fixed << std::setprecision(6) << "bra imag:\n" << snxbraIA.imag() << "\n";
			std::cout << "\n";
			auto snxketIA = readMatrix("k" + std::to_string(I) + cartDict[alpha]);
			std::cout << std::fixed << std::setprecision(6) << "ket real:\n" << snxketIA.real() << "\n";
			std::cout << std::fixed << std::setprecision(6) << "ket imag:\n" << snxketIA.imag() << "\n";

			


			//std::cout << std::fixed << std::setprecision(6) << "real:\n" << snxbraSum.real() << "\n";
			//std::cout << std::fixed << std::setprecision(6) << "imag:\n" << snxbraSum.imag() << "\n";
			//std::cout << std::fixed << std::setprecision(6) << "real:\n" << snxketSum.real() << "\n";
			//std::cout << std::fixed << std::setprecision(6) << "imag:\n" << snxketSum.imag() << "\n";

			std::cout << "\n\n\n";
		}
	}
	//*/





	// calculate and decompose density matrix
	//Eigen::MatrixXcd denMat = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	//for (int k=0; k<spinorSize; k++) {
	//	for (int l=0; l<spinorSize; l++) {
	//		for (int i=0; i<nocc; i++) {
	//			denMat(k, l) += std::conj(spinor(k, i)) * spinor(l, i);
	//		}
	//	}
	//}
	//std::cout << "\ndensity matrix:\n" << denMat << "\n\n";
	//const auto aa = denMat(Eigen::seq(0,spinorSize/2-1), Eigen::seq(0,spinorSize/2-1));
	//const auto ab = denMat(Eigen::seq(0,spinorSize/2-1), Eigen::seq(spinorSize/2, spinorSize-1));
	//const auto ba = denMat(Eigen::seq(spinorSize/2, spinorSize-1), Eigen::seq(0,spinorSize/2-1));
	//const auto bb = denMat(Eigen::seq(spinorSize/2, spinorSize-1), Eigen::seq(spinorSize/2, spinorSize-1));
	//std::cout << "\nalpha-alpha block:\n" << aa << "\n\n";
	//std::cout << "\nbeta-beta block:\n" << bb << "\n\n";
	//std::cout << "\nalpha-beta block:\n" << ab << "\n\n";
	//std::cout << "\nbeta-alpha block:\n" << ba << "\n\n";
	//std::cout << "\n0-Komponente:\n" << (aa+bb)/2 << "\n\n";
	//std::cout << "\n1-Komponente (x):\n" << (ab+ba)/2 << "\n\n";
	//std::cout << "\n2-Komponente (y):\n" << 1.0i*(ab-ba)/2 << "\n\n";
	//std::cout << "\n3-Komponente (z):\n" << (aa-bb)/2 << "\n\n";

	
	

	/* debug derivatives of coefficients
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			std::cout << "\nperturbed coefficients for atom " << I << cartDict[alpha] << ":\n";
			Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
			const auto snxbra = readHerm("b" + std::to_string(I) + cartDict[alpha]);
			const auto snxket = readHerm("k" + std::to_string(I) + cartDict[alpha]);
			const auto snx = snxbra + snxket;
			Eigen::MatrixXcd snxBig(spinorSize, spinorSize);
			snxBig << snx, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
				Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snx;
			const auto snxMO = spinor.adjoint() * snxBig * spinor;

			for (int i=0; i<nocc; i++) {
				for (int j=0; j<nocc; j++) {
					U(i, j) = -0.5 * snxMO(i, j);
				}
			}
			const auto uIA = readVector("u" + std::to_string(I) + "_" + std::to_string(alpha));
			for (int i=0; i<nocc; i++) {
				for (int a=nocc; a<spinorSize; a++) {
					U(a, i) = uIA(i*nvirt+a-nocc);
					U(i, a) = std::conj(uIA(i*nvirt+a-nocc));
				}
			}

			Eigen::MatrixXcd cnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
			for (int mu=0; mu<spinorSize; mu++) {
				for (int i=0; i<spinorSize; i++) {
					for (int q=0; q<spinorSize; q++) {
						cnx(mu, i) += spinor(mu, q) * U(q, i);
					}
				}
			}
			//std::cout << std::setprecision(4) << cnx << "\n\n";
			
		}
	}

	//std::cout << "\nunperturbed spinor:\n" << spinor << "\n\n";
	//*/






	/* =============================================================================================== 1
	{int I = 0;
	int alpha = 0;
	int J = 0;
	int beta = 0;
	std::cout << "\nperturbed coefficients for atom " << I << cartDict[alpha] << ":\n";
	Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	auto snxbra = readHerm("b" + std::to_string(I) + cartDict[alpha]);
	auto snxket = readHerm("k" + std::to_string(I) + cartDict[alpha]);
	auto snx = snxbra + snxket;
	Eigen::MatrixXcd snxBig(spinorSize, spinorSize);
	snxBig << snx, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snx;
	const auto snxMO = spinor.adjoint() * snxBig * spinor;

	// switch lower triangle of snxbra and snxket
	auto tmpbra = snxbra;
	for (int i=0; i<spinorSize/2; i++) {
		for (int j=0; j<i; j++) {
			snxbra(i, j) = snxket(i, j);
			snxket(i, j) = tmpbra(i, j);
		}
	}

	for (int i=0; i<nocc; i++) {
		for (int j=0; j<nocc; j++) {
			U(i, j) = -0.5 * snxMO(i, j);
		}
	}
	const auto uIA = readVector("u" + std::to_string(I) + "_" + std::to_string(alpha));
	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			U(a, i) = uIA(i*nvirt+a-nocc);
			U(i, a) = std::conj(uIA(i*nvirt+a-nocc));
		}
	}

	Eigen::MatrixXcd cnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	for (int mu=0; mu<spinorSize; mu++) {
		for (int i=0; i<spinorSize; i++) {
			for (int q=0; q<spinorSize; q++) {
				cnx(mu, i) += spinor(mu, q) * U(q, i);
			}
		}
	}


	//std::cout << "\nperturbed coefficients for atom " << J << cartDict[beta] << ":\n";
	Eigen::MatrixXcd U2 = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	auto snxbra2 = readHerm("b" + std::to_string(J) + cartDict[beta]);
	auto snxket2 = readHerm("k" + std::to_string(J) + cartDict[beta]);
	auto snx2 = snxbra2 + snxket2;
	Eigen::MatrixXcd snxBig2(spinorSize, spinorSize);
	snxBig2 << snx2, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snx2;
	const auto snxMO2 = spinor.adjoint() * snxBig2 * spinor;

	// switch lower triangle of snxbra2 and snxket2
	auto tmpbra2 = snxbra2;
	for (int i=0; i<spinorSize/2; i++) {
		for (int j=0; j<i; j++) {
			snxbra2(i, j) = snxket2(i, j);
			snxket2(i, j) = tmpbra2(i, j);
		}
	}

	for (int i=0; i<nocc; i++) {
		for (int j=0; j<nocc; j++) {
			U2(i, j) = -0.5 * snxMO2(i, j);
		}
	}
	const auto uIA2 = readVector("u" + std::to_string(J) + "_" + std::to_string(beta));
	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			U2(a, i) = uIA2(i*nvirt+a-nocc);
			U2(i, a) = std::conj(uIA2(i*nvirt+a-nocc));
		}
	}

	Eigen::MatrixXcd cnx2 = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	for (int mu=0; mu<spinorSize; mu++) {
		for (int i=0; i<spinorSize; i++) {
			for (int q=0; q<spinorSize; q++) {
				cnx2(mu, i) += spinor(mu, q) * U2(q, i);
			}
		}
	}


	// overlap stuff
	auto smat = readHerm("smat");
	Eigen::MatrixXcd smatBig = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	smatBig << smat, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), smat;


	Eigen::MatrixXcd smatBra(spinorSize, spinorSize);
	smatBra << snxbra, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbra;

	Eigen::MatrixXcd smatKet(spinorSize, spinorSize);
	smatKet << snxket2, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxket2;

	const auto braket = readHerm("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta]);
	Eigen::MatrixXcd tmp3(spinorSize, spinorSize);
	tmp3 << braket, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braket;


	//auto kek = cnx.adjoint() * smatBig * cnx2 + cnx.adjoint() * smatKet * spinor + spinor.adjoint() * smatBra * cnx2 + spinor.adjoint() * tmp3 * spinor;

	//std::cout << "\nlkasjdfölaksdjf:\n" << smatBig << "\n\n";

	std::complex<double> res = (0, 0);
	for (int i=0; i<nocc; i++) {
		//res += (cnx.adjoint() * smatBig * cnx2)(i, i);
		//res += (cnx.adjoint()    * smatKet * spinor)(i, i);
		//res += (spinor.adjoint() * smatBra * cnx2)(i, i);
		//res += (spinor.adjoint() * tmp3 * spinor)(i, i);
	}

	//std::cout << "\nerster Beitrag:\n" << smatKet << "\n\n";
	//std::cout << "\nzweiter Beitrag:\n" << smatBra << "\n\n";
	//std::cout << "\nerster Beitrag:\n"  << (cnx.adjoint()    * smatKet * spinor) << "\n\n";
	//std::cout << "\nzweiter Beitrag:\n" << (spinor.adjoint() * smatBra * cnx2)   << "\n\n";

	for (int j=0; j<nocc; j++) {
		for (int mu=0; mu<spinorSize; mu++) {
			for (int nu=0; nu<spinorSize; nu++) {
				res += std::conj(cnx(mu, j)) * smatBig(mu, nu) * cnx2(nu, j);
				res += std::conj(cnx(mu, j)) * smatKet(mu, nu) * spinor(nu, j);
				res += std::conj(spinor(mu, j)) * smatBra(mu, nu) * cnx2(nu, j);
				res += std::conj(spinor(mu, j)) * tmp3(mu, nu) * spinor(nu, j);
			}
		}
	}

	std::cout << " berry? = " << res << "\n";}
	//*/


	/**************************************************************************************************************
	 *                                             Debug von cnx                                                  *
	 *************************************************************************************************************/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///*
	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			std::cout << "\n\n\n\nDebug atom " << I << ", cart " << alpha << "\n\n";
			//if (argc>1) std::cout << "\n\n\n\nDebug atom " << I << ", cart " << alpha << "\n\n";
			//Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
			//auto snxbra = readHerm("b" + std::to_string(I) + cartDict[alpha]);
			//auto snxket = readHerm("k" + std::to_string(I) + cartDict[alpha]);
			//auto snx = snxbra + snxket;
			//Eigen::MatrixXcd snxBig(spinorSize, spinorSize);
			//snxBig << snx, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
			//	Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snx;
			//const auto snxMO = spinor.adjoint() * snxBig * spinor;

			//std::cout << std::fixed << std::setprecision(8) << "\nsnxbra VOR SWAP:\n" << snxbra << "\n\n";
			//std::cout << std::fixed << std::setprecision(8) << "\nsnxket VOR SWAP:\n" << snxket << "\n\n";

			//// switch lower triangle of snxbra and snxket
			//auto tmpbra = snxbra;
			//for (int i=0; i<spinorSize/2; i++) {
			//	for (int j=0; j<i; j++) {
			//		snxbra(i, j) = snxket(i, j);
			//		snxket(i, j) = tmpbra(i, j);
			//	}
			//}

			//std::cout << std::fixed << std::setprecision(8) << "\nsnxbra NACH SWAP:\n" << snxbra << "\n\n";
			//std::cout << std::fixed << std::setprecision(8) << "\nsnxket NACH SWAP:\n" << snxket << "\n\n";

			//Eigen::MatrixXcd snxBraBig(spinorSize, spinorSize);
			//snxBraBig << snxbra, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
			//	Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbra;
			//std::cout << std::fixed << std::setprecision(8) << "\nsnxbraBig:\n" << snxBraBig << "\n\n";
			//const auto snxBraMO = spinor.adjoint() * snxBraBig * spinor;
			//std::cout << std::fixed << std::setprecision(8) << "\nsnxbraMO:\n" << snxBraMO << "\n\n";
			//if (argc>1) { if (*argv[1]=='0') { std::cout << std::fixed << std::setprecision(8) << "\nsnxbraMO:\n" << snxBraMO << "\n\n"; } }
			//std::cout << std::fixed << std::setprecision(8) << "\nsnxbra:\n" << snxbra << "\n\n";

			//Eigen::MatrixXcd snxKetBig(spinorSize, spinorSize);
			//snxKetBig << snxket, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
			//	Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxket;
			//const auto snxKetMO = spinor.adjoint() * snxKetBig * spinor;
			//std::cout << std::fixed << std::setprecision(8) << "\nsnxketMO:\n" << snxKetMO << "\n\n";
			//if (argc>1) { if (*argv[1]=='1') { std::cout << std::fixed << std::setprecision(8) << "\nsnxketMO:\n" << snxKetMO << "\n\n"; } }
			//std::cout << std::fixed << std::setprecision(8) << "\nsnxket:\n" << snxket << "\n\n";

			//for (int i=0; i<nocc; i++) {
			//	for (int j=0; j<nocc; j++) {
			//		U(i, j) = -0.5 * snxMO(i, j);
			//	}
			//}
			const auto uIA = readVector("u" + std::to_string(I) + "_" + std::to_string(alpha));
			//for (int i=0; i<nocc; i++) {
			//	for (int a=nocc; a<spinorSize; a++) {
			//		U(a, i) = uIA(i*nvirt+a-nocc);
			//		U(i, a) = std::conj( -1.0*(snxMO(a, i) + uIA(i*nvirt+a-nocc)) );
			//	}
			//}


			//Eigen::MatrixXcd cnxAna = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
			//for (int mu=0; mu<spinorSize; mu++) {
			//	for (int i=0; i<spinorSize; i++) {
			//		for (int q=0; q<spinorSize; q++) {
			//			cnxAna(mu, i) += spinor(mu, q) * U(q, i);
			//		}
			//	}
			//}

			auto cplus  = readNumSpinor("TEST5/" + std::to_string(3*I + alpha+1) + "IPlus.out");
			auto cminus = readNumSpinor("TEST5/" + std::to_string(3*I + alpha+1) + "IMinus.out");
			auto cnxNum = (cplus - cminus) / 2e-3;
			auto uNum = spinor.inverse() * cnxNum;
			//std::cout << std::fixed << std::setprecision(8) << uNum << "\n";
			//std::cout << "\n\n";
			//std::cout << std::fixed << std::setprecision(8) << U << "\n";

			// calc bNum
			Eigen::VectorXcd uIANum(2*nocc*nvirt);
			for (int i=0; i<nocc; i++) {
				for (int a=nocc; a<spinorSize; a++) {
					//U(a, i) = uIA(i*nvirt+a-nocc);
					uIANum(i*nvirt+a-nocc) = uNum(a, i);
					uIANum(nocc*nvirt + i*nvirt+a-nocc) = std::conj(uNum(a, i));
				}
			}
			Eigen::MatrixXcd H(2*nocc*nvirt, 2*nocc*nvirt);
			//H << A, B, B.conjugate(), A.conjugate();
			H << A.conjugate(), B.conjugate(), B, A;
			const auto bNum = H*uIANum;
			const auto b = readVector("b0ai" + std::to_string(I) + "_" + std::to_string(alpha));

			//std::cout << "snxBraMo * UNum:\n" << snxBraMO*uNum << "\n\n";

			//std::cout << std::fixed << std::setprecision(8) << "bNum:\n" << bNum << "\n";
			//std::cout << "\n\n";
			//std::cout << std::fixed << std::setprecision(8) << "b:\n" << b << "\n";
			std::cout << "      uNum			uAnalytic\n";
			for (int i=0; i<nocc; i++) for (int a=nocc; a<spinorSize; a++) std::cout << std::fixed << std::setprecision(10) << uNum(a, i) << "		" << uIA(i*nvirt + a - nocc) << "\n";
			std::cout << "\n\n";
			std::cout << "      bNum			bAnalytic\n";
			for (int i=0; i<nocc*nvirt; i++) {std::cout << std::fixed << std::setprecision(10) << bNum(i) << "		" << b(i) << "\n";}
			std::cout << "\n\n";

			//std::cout << "Num                  Ana\n";
			//for (int i=0; i<spinorSize; i++) {
			//	for (int j=0; j<spinorSize; j++) {
			//		std::cout << std::fixed << std::setprecision(8) << cnxNum(i, j) << "			" << cnxAna(i, j) << "\n";
			//	}
			//}
			//std::cout << std::fixed << std::setprecision(8) << "num:\n" << cnxNum << "\n";
			//std::cout << "\n\n";
			//std::cout << std::fixed << std::setprecision(8) << "ana:\n" << cnxAna << "\n";
			//std::cout << "\n\n\n\n";
			
			// fock debug
			//auto cplus  = readNumSpinor("TEST5/" + std::to_string(alpha+1) + "IPlus.out");
			//auto cminus = readNumSpinor("TEST5/" + std::to_string(alpha+1) + "IMinus.out");
			//const auto hmat = readHerm("hmat");
			//const auto hmatplus = readHerm(std::to_string(I) + "_" + cartDict[alpha] + "plus/hmat");
		}
	}
	//*/
	



	/*
	std::cout << "\n\n\n   alternativer Weg...\n\n\n";
	Eigen::MatrixXcd berry2a(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry2b(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry2b2(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry2c(3*atomNum, 3*atomNum);
	Eigen::MatrixXcd berry2d(3*atomNum, 3*atomNum);

	for (int I=0; I<atomNum; I++) {
		for (int alpha=0; alpha<3; alpha++) {
			for (int J=0; J<atomNum; J++) {
				for (int beta=0; beta<3; beta++) 
	// ================================================================================================ 2
	{
	//int I = 0;
	//int alpha = 0;
	//int J = 0;
	//int beta = 0;
	Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	auto snxbra = readHerm("b" + std::to_string(I) + cartDict[alpha]);
	auto snxket = readHerm("k" + std::to_string(I) + cartDict[alpha]);
	auto snx = snxbra + snxket;
	Eigen::MatrixXcd snxBig(spinorSize, spinorSize);
	snxBig << snx, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snx;
	const auto snxMO = spinor.adjoint() * snxBig * spinor;

	// switch lower triangle of snxbra and snxket
	auto tmpbra = snxbra;
	for (int i=0; i<spinorSize/2; i++) {
		for (int j=0; j<i; j++) {
			snxbra(i, j) = snxket(i, j);
			snxket(i, j) = tmpbra(i, j);
		}
	}

	for (int i=0; i<nocc; i++) {
		for (int j=0; j<nocc; j++) {
			U(i, j) = -0.5 * snxMO(i, j);
		}
	}
	const auto uIA = readVector("u" + std::to_string(I) + "_" + std::to_string(alpha));
	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			U(a, i) = uIA(i*nvirt+a-nocc);
			U(i, a) = std::conj( -1.0*(snxMO(a, i) + uIA(i*nvirt+a-nocc)) );
			//U(i, a) = uIA(i*nvirt+a-nocc);
			//U(a, i) = std::conj( -1.0*(snxMO(i, a) + uIA(i*nvirt+a-nocc)) );
		}
	}

	//Eigen::MatrixXcd cnx = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	//for (int mu=0; mu<spinorSize; mu++) {
	//	for (int i=0; i<spinorSize; i++) {
	//		for (int q=0; q<spinorSize; q++) {
	//			cnx(mu, i) += spinor(mu, q) * U(q, i);
	//		}
	//	}
	//}
	auto cplus  = readNumSpinor("TEST5/" + std::to_string(3*I + alpha+1) + "IPlus.out");
	auto cminus = readNumSpinor("TEST5/" + std::to_string(3*I + alpha+1) + "IMinus.out");
	auto cnx = (cplus - cminus) / 2e-3;

	//std::cout << "\nperturbed coefficients for atom " << I << cartDict[alpha] << ":\n" << spinor + 1e-3*cnx << "\n\n";
	//std::cout << "\nU for atom " << I << cartDict[alpha] << ":\n" << U << "\n\n";
	//std::cout << "\nU * U^H for atom " << I << cartDict[alpha] << ":\n" << U*U.adjoint() << "\n\n";


	//std::cout << "\nperturbed coefficients for atom " << J << cartDict[beta] << ":\n";
	Eigen::MatrixXcd U2 = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	auto snxbra2 = readHerm("b" + std::to_string(J) + cartDict[beta]);
	auto snxket2 = readHerm("k" + std::to_string(J) + cartDict[beta]);
	auto snx2 = snxbra2 + snxket2;
	Eigen::MatrixXcd snxBig2(spinorSize, spinorSize);
	snxBig2 << snx2, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snx2;
	const auto snxMO2 = spinor.adjoint() * snxBig2 * spinor;

	//std::cout << "\nsnxbra2:\n" << snxbra2 << "\n\n";
	//std::cout << "\nsnxket2:\n" << snxket2 << "\n\n";

	// switch lower triangle of snxbra2 and snxket2
	auto tmpbra2 = snxbra2;
	for (int i=0; i<spinorSize/2; i++) {
		for (int j=0; j<i; j++) {
			snxbra2(i, j) = snxket2(i, j);
			snxket2(i, j) = tmpbra2(i, j);
		}
	}
	//std::cout << "\nsnxbra2:\n" << snxbra2 << "\n\n";
	//std::cout << "\nsnxket2:\n" << snxket2 << "\n\n";

	for (int i=0; i<nocc; i++) {
		for (int j=0; j<nocc; j++) {
			U2(i, j) = -0.5 * snxMO2(i, j);
		}
	}
	const auto uIA2 = readVector("u" + std::to_string(J) + "_" + std::to_string(beta));
	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			U2(a, i) = uIA2(i*nvirt+a-nocc);
			U2(i, a) = std::conj( -1.0*(snxMO2(a, i) + uIA2(i*nvirt+a-nocc)) );
			//U2(i, a) = uIA2(i*nvirt+a-nocc);
			//U2(a, i) = std::conj( -1.0*(snxMO2(i, a) + uIA2(i*nvirt+a-nocc)) );
		}
	}

	//Eigen::MatrixXcd cnx2 = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	//for (int mu=0; mu<spinorSize; mu++) {
	//	for (int i=0; i<spinorSize; i++) {
	//		for (int q=0; q<spinorSize; q++) {
	//			cnx2(mu, i) += spinor(mu, q) * U2(q, i);
	//		}
	//	}
	//}
	auto cplus2  = readNumSpinor("TEST5/" + std::to_string(3*J + beta+1) + "IPlus.out");
	auto cminus2 = readNumSpinor("TEST5/" + std::to_string(3*J + beta+1) + "IMinus.out");
	auto cnx2 = (cplus2 - cminus2) / 2e-3;



	// overlap stuff
	auto smat = readHerm("smat");
	Eigen::MatrixXcd smatBig = Eigen::MatrixXcd::Zero(spinorSize, spinorSize);
	smatBig << smat, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), smat;


	Eigen::MatrixXcd smatBra(spinorSize, spinorSize);
	smatBra << snxbra, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxbra;

	Eigen::MatrixXcd smatKet(spinorSize, spinorSize);
	smatKet << snxket2, Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), snxket2;

	// doppelt abgeleitete overlap matrix
	auto braketA = readHerm("bk" + std::to_string(I) + cartDict[alpha] + std::to_string(J) + cartDict[beta]);
	auto braketB = readHerm("bk" + std::to_string(J) + cartDict[beta] + std::to_string(I) + cartDict[alpha]);
	// switch lower triangle
	auto tmpbraket = braketA;
	for (int i=0; i<spinorSize/2; i++) {
		for (int j=0; j<i; j++) {
			braketA(i, j) = braketB(i, j);
			braketB(i, j) = tmpbraket(i, j);
		}
	}
	Eigen::MatrixXcd tmp3(spinorSize, spinorSize);
	tmp3 << braketA.conjugate(), Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2),
		Eigen::MatrixXcd::Zero(spinorSize/2, spinorSize/2), braketA.conjugate();


	//auto kek = cnx.adjoint() * smatBig * cnx2 + cnx.adjoint() * smatKet * spinor + spinor.adjoint() * smatBra * cnx2 + spinor.adjoint() * tmp3 * spinor;

	std::complex<double> resa = (0, 0);
	std::complex<double> resb = (0, 0);
	std::complex<double> resb2 = (0, 0);
	std::complex<double> resc = (0, 0);
	for (int i=0; i<nocc; i++) {
		resa += (cnx.adjoint() * smatBig * cnx2)(i, i);
		resb += (cnx.adjoint()    * smatKet * spinor)(i, i);
		resb2 += (spinor.adjoint() * smatBra * cnx2)(i, i);
		resc += (spinor.adjoint() * tmp3 * spinor)(i, i);
	}

	//std::cout << "\nerster Beitrag:\n"  << smatKet << "\n\n";
	//std::cout << "\nzweiter Beitrag:\n" << smatBra   << "\n\n";

	for (int j=0; j<nocc; j++) {
		for (int mu=0; mu<spinorSize; mu++) {
			for (int nu=0; nu<spinorSize; nu++) {
				//res += std::conj(cnx(mu, j)) * smatBig(mu, nu) * cnx2(nu, j);
				//res += std::conj(cnx(mu, j)) * smatKet(mu, nu) * spinor(nu, j);
				//res += std::conj(spinor(mu, j)) * smatBra(mu, nu) * cnx2(nu, j);
				//res += std::conj(spinor(mu, j)) * tmp3(mu, nu) * spinor(nu, j);
			}
		}
	}

	//std::cout << "\ncnx:\n" << cnx << "\n\n";
	//std::cout << "\ncnx2:\n" << cnx2 << "\n\n";

	//std::cout << "\nsnxbra:\n" << snxbra << "\n\n";
	//std::cout << "\nsnxbra2:\n" << snxbra2 << "\n\n";
	//std::cout << "\nsnxket:\n" << snxket << "\n\n";
	//std::cout << "\nsnxket2:\n" << snxket2 << "\n\n";

	berry2a(3*I+alpha, 3*J+beta) = -2.0*resa;
	berry2b(3*I+alpha, 3*J+beta) = -2.0*resb;
	berry2b2(3*I+alpha, 3*J+beta) = -2.0*resb2;
	berry2c(3*I+alpha, 3*J+beta) = -2.0*resc;
	berry2d(3*I+alpha, 3*J+beta) = -2.0*resa - 2.0*resb - 2.0*resb2 - 2.0*resc;
	//std::cout << " berry? = " << res << "\n";
	}}}}
	//std::cout << "\nberry2:\n" << berry2 << "\n\n";
	//std::cout << "\n\n\nBerry-curvature sym:\n" << 0.5*(berry2 + berry2.transpose()) << "\n\n";
	//std::cout << "\n\n\nBerry-curvature antisym:\n" << 0.5*(berry2 - berry2.transpose()) << "\n\n";
	std::cout << std::fixed;
	std::cout << std::setprecision(10);
	std::cout << "\n\n" << berry2d.imag() << "\n\n";
	//std::cout << "\n\nantisym imag part C' S C':\n" << 0.5*(berry2a - berry2a.transpose()).imag() << "\n\n";
	//std::cout << "\n\nantisym imag part C' S' C + C S' C':\n" << 0.5*(berry2b - berry2b.transpose()).imag() << "\n\n";
	//std::cout << "\n\nantisym imag part C' S' C:\n" << 0.5*(berry2b - berry2b.transpose()).imag() << "\n\n";
	//std::cout << "\n\nantisym imag part C S' C':\n" << 0.5*(berry2b2 - berry2b2.transpose()).imag() << "\n\n";
	//std::cout << "\n\nantisym imag part C S'' C:\n" << 0.5*(berry2c - berry2c.transpose()).imag() << "\n\n";
	//std::cout << "\n\ntotal:\n" << 0.5*(berry2d - berry2d.transpose()).imag() << "\n\n";
	//*/
	

	return 0;
}
