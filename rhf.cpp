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

int main() {
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
	const auto mos = readMos(epsilon);
	std::cout << "\nmos:\n" << mos << "\n\n";

	const int nAO = mos.cols();

	Eigen::MatrixXcd dmat = Eigen::MatrixXcd::Zero(nAO, nAO);
	for (int mu=0; mu<nAO; mu++) {
		for (int nu=0; nu<nAO; nu++) {
			for (int i=0; i<nocc/2; i++) {
				dmat(mu, nu) += 2.0 * mos(mu, i) * mos(nu, i);
			}
		}
	}
	std::cout << "\ndmat:\n" << dmat << "\n\n";

	/**************************************************************************************************************
	 *                                               fock matrix                                                  *
	 *************************************************************************************************************/
	const auto smat = readHerm("smat");
	std::cout << "\nsmat:\n" << smat << "\n\n";
	const auto smatMO = mos.transpose() * smat * mos;
	std::cout << "\nsmatMO:\n" << smatMO << "\n\n";

	const auto hmat = readHerm("hmat");
	std::cout << "\nhmat:\n" << hmat << "\n\n";
	const auto hmatMO = mos.transpose() * hmat * mos;
	std::cout << "\nhmatMO:\n" << hmatMO << "\n\n";

	auto fci = readFourcenter("");
	Eigen::MatrixXcd jmat = Eigen::MatrixXcd::Zero(nAO, nAO);
	Eigen::MatrixXcd xmat = Eigen::MatrixXcd::Zero(nAO, nAO);

	for (int mu=0; mu<nAO; mu++) {
		for (int nu=0; nu<nAO; nu++) {
			for (int ka=0; ka<nAO; ka++) {
				for (int la=0; la<nAO; la++) {
					jmat(mu, nu) += 1.0 * dmat(ka, la) * fci[mu][nu][ka][la];
					xmat(mu, nu) -= 0.5 * dmat(ka, la) * fci[mu][ka][nu][la];
				}
			}
		}
	}
	std::cout << "\njmat:\n" << jmat << "\n\n";
	const auto jmatMO = mos.transpose() * jmat * mos;
	std::cout << "\njmatMO:\n" << jmatMO << "\n\n";
	std::cout << "\nxmat:\n" << xmat << "\n\n";
	const auto xmatMO = mos.transpose() * xmat * mos;
	std::cout << "\nxmatMO:\n" << xmatMO << "\n\n";

	const auto fmat = hmat + jmat + xmat;
	std::cout << "\nfmat:\n" << fmat << "\n\n";
	const auto fmatMO = mos.transpose() * fmat * mos;
	std::cout << "\nfmatMO:\n" << fmatMO << "\n\n";



	/**************************************************************************************************************
	 *                                            d fock / d Nx (Num)                                             *
	 *************************************************************************************************************/
	std::cout << "\n\n\n\nfnx stuff\n\n";

	const auto smatplus  = readHerm("x0plus/smat");
	const auto smatminus = readHerm("x0minus/smat");
	const auto snxNum = (smatplus - smatminus) / 2e-4;
	std::cout << "\nsnxNum (SAO):\n" << snxNum << "\n\n";
	const auto snxNumMO = mos.transpose() * snxNum * mos;
	std::cout << "\nsnxNum (MO):\n" << snxNumMO << "\n\n";

	const auto hmatplus  = readHerm("x0plus/hmat");
	const auto hmatminus = readHerm("x0minus/hmat");
	const auto hnxNum = (hmatplus - hmatminus) / 2e-4;
	std::cout << "\nhnxNum (SAO):\n" << hnxNum << "\n\n";
	const auto hnxNumMO = mos.transpose() * hnxNum * mos;
	std::cout << "\nhnxNum (MO):\n" << hnxNumMO << "\n\n";

	auto fciplus  = readFourcenter("x0plus/");
	auto fciminus = readFourcenter("x0minus/");
	fourD fcinxNum(nAO,
			std::vector<std::vector<std::vector<std::complex<double>>>>(nAO,
				std::vector<std::vector<std::complex<double>>>(nAO,
					std::vector<std::complex<double>>(nAO))));
	for (int i=0; i<nAO; i++) {
		for (int j=0; j<nAO; j++) {
			for (int k=0; k<nAO; k++) {
				for (int l=0; l<nAO; l++) {
					fcinxNum[i][j][k][l] = (fciplus[i][j][k][l] - fciminus[i][j][k][l]) / 2e-4;
				}
			}
		}
	}

	Eigen::MatrixXcd jnxNum = Eigen::MatrixXcd::Zero(nAO, nAO);
	Eigen::MatrixXcd xnxNum = Eigen::MatrixXcd::Zero(nAO, nAO);

	for (int mu=0; mu<nAO; mu++) {
		for (int nu=0; nu<nAO; nu++) {
			for (int ka=0; ka<nAO; ka++) {
				for (int la=0; la<nAO; la++) {
					jnxNum(mu, nu) += 1.0 * dmat(ka, la) * fcinxNum[mu][nu][ka][la];
					xnxNum(mu, nu) -= 0.5 * dmat(ka, la) * fcinxNum[mu][ka][nu][la];
				}
			}
		}
	}
	std::cout << "\njnxNum:\n" << jnxNum << "\n\n";
	std::cout << "\nxnxNum:\n" << xnxNum << "\n\n";

	const auto fnxNum = hnxNum + jnxNum + xnxNum;
	std::cout << "\nfnxNum:\n" << fnxNum << "\n\n";
	const auto fnxNumMO = mos.transpose() * fnxNum * mos;
	std::cout << "\nfnxNumMO:\n" << fnxNumMO << "\n\n";


	/**************************************************************************************************************
	 *                                              fci to MO basis                                               *
	 *************************************************************************************************************/
	std::cout << "\n:: transforming fourcenter integrals to MO basis...\n" << std::flush;
	fourD fourCenterIntegral(nAO,
			std::vector<std::vector<std::vector<std::complex<double>>>>(nAO,
				std::vector<std::vector<std::complex<double>>>(nAO,
					std::vector<std::complex<double>>(nAO))));

	fourD tmp(nAO,
			std::vector<std::vector<std::vector<std::complex<double>>>>(nAO,
				std::vector<std::vector<std::complex<double>>>(nAO,
					std::vector<std::complex<double>>(nAO))));
	fourD tmp2(nAO,
			std::vector<std::vector<std::vector<std::complex<double>>>>(nAO,
				std::vector<std::vector<std::complex<double>>>(nAO,
					std::vector<std::complex<double>>(nAO))));
	fourD tmp3(nAO,
			std::vector<std::vector<std::vector<std::complex<double>>>>(nAO,
				std::vector<std::vector<std::complex<double>>>(nAO,
					std::vector<std::complex<double>>(nAO))));


	std::cout << "      first transformation..." << std::flush;
	// first transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<nAO; i++) {
		for (int nu=0; nu<nAO; nu++) {
			for (int lam=0; lam<nAO; lam++) {
				for (int sig=0; sig<nAO; sig++) {
					tmp[i][nu][lam][sig] = 0;
					for (int mu=0; mu<nAO; mu++) {
						tmp[i][nu][lam][sig] += mos(mu, i) * fci[mu][nu][lam][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;

	std::cout << "      second transformation..." << std::flush;
	// second transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<nAO; i++) {
		for (int j=0; j<nAO; j++) {
			for (int lam=0; lam<nAO; lam++) {
				for (int sig=0; sig<nAO; sig++) {
					tmp2[i][j][lam][sig] = 0;
					for (int nu=0; nu<nAO; nu++) {
						tmp2[i][j][lam][sig] += mos(nu, j) * tmp[i][nu][lam][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;

	std::cout << "      third transformation..." << std::flush;
	// third transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<nAO; i++) {
		for (int j=0; j<nAO; j++) {
			for (int k=0; k<nAO; k++) {
				for (int sig=0; sig<nAO; sig++) {
					tmp3[i][j][k][sig] = 0;
					for (int lam=0; lam<nAO; lam++) {
						tmp3[i][j][k][sig] += mos(lam, k) * tmp2[i][j][lam][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;

	std::cout << "      fourth transformation..." << std::flush;
	// fourth transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<nAO; i++) {
		for (int j=0; j<nAO; j++) {
			for (int k=0; k<nAO; k++) {
				for (int l=0; l<nAO; l++) {
					fourCenterIntegral[i][j][k][l] = 0;
					for (int sig=0; sig<nAO; sig++) {
						fourCenterIntegral[i][j][k][l] += mos(sig, l) * tmp3[i][j][k][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;
	std::cout << "\n" << std::flush;



	/**************************************************************************************************************
	 *                                                     b0ai                                                   *
	 *************************************************************************************************************/
	Eigen::VectorXcd b0ai = Eigen::VectorXcd::Zero(nocc*nvirt/4);
	std::cout << "\n\n\n\nb0ai nur mit 1. Term:\n";
	for (int i=0; i<nocc/2; i++) {
		for (int a=nocc/2; a<nAO; a++) {
			const int index = i*(nvirt/2) + a - (nocc/2);
			b0ai(index) = -fnxNumMO(a, i);
		}
	}
	std::cout << "\n" << b0ai << "\n\n";

	std::cout << "\n\n\n\nb0ai nur mit 1. + 2. Term:\n";
	for (int i=0; i<nocc/2; i++) {
		for (int a=nocc/2; a<nAO; a++) {
			const int index = i*(nvirt/2) + a - (nocc/2);
			for (int j=0; j<nocc/2; j++) {
				b0ai(index) += snxNumMO(a, j) * fmatMO(j, i);
			}
		}
	}
	std::cout << "\n" << b0ai << "\n\n";
	
	std::cout << "\n\n\n\nb0ai nur mit allen Termen:\n";
	for (int i=0; i<nocc/2; i++) {
		for (int a=nocc/2; a<nAO; a++) {
			const int index = i*(nvirt/2) + a - (nocc/2);
			for (int k=0; k<nocc/2; k++) {
				for (int l=0; l<nocc/2; l++) {
					b0ai(index) += 2.0 * snxNumMO(k, l) * fourCenterIntegral[a][i][l][k];
					b0ai(index) -= snxNumMO(k, l) * fourCenterIntegral[a][k][l][i];
				}
			}
		}
	}
	std::cout << "\n" << b0ai << "\n\n";
	
}
