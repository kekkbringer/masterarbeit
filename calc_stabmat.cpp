#include "calc_stabmat.hpp"
#include "read_fourcenter.hpp"
#include "read_spinor.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <iomanip>

void calcStabmat(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B) {
    // looking for occupied spinors in control file
	std::cout << "\n:: reading control file...\n" << std::flush;
	int nocc = 0;
	std::ifstream control("control");
	if (control.fail()) {std::cout << "\nWARNING: cant read control file!\n";}
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
		}
    }


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


    // reading four center integrals... oh boi...
	std::cout << " reading four center integrals...\n" << std::flush;
	const auto fci = readFourcenter();
	std::cout << "      size of fci:   " << fci.size() << " x " << fci.size() << " x " << fci.size() << " x " << fci.size() << "\n";
	std::cout << "      amounts to:    " << (sizeof fci)*fci.size()*fci.size()*fci.size()*fci.size()/1000000.0 << " MB\n";
	std::cout << "\n" << std::flush;




    // fourcenter integral transformation into MO-basis
	std::cout << "\n:: transforming fourcenter integrals to spinor basis...\n" << std::flush;
	fourD fourCenterIntegral(spinorSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
				std::vector<std::vector<std::complex<double>>>(spinorSize,
					std::vector<std::complex<double>>(spinorSize))));

	fourD tmp(spinorSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
				std::vector<std::vector<std::complex<double>>>(spinorSize,
					std::vector<std::complex<double>>(spinorSize))));
	fourD tmp2(spinorSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
				std::vector<std::vector<std::complex<double>>>(spinorSize,
					std::vector<std::complex<double>>(spinorSize))));
	fourD tmp3(spinorSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
				std::vector<std::vector<std::complex<double>>>(spinorSize,
					std::vector<std::complex<double>>(spinorSize))));

	// double dim of fci
	fourD fciD(spinorSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(spinorSize,
				std::vector<std::vector<std::complex<double>>>(spinorSize,
					std::vector<std::complex<double>>(spinorSize))));

	std::cout << "      size of new fci:   " << fciD.size() << " x " << fciD.size() << " x " << fciD.size() << " x " << fciD.size() << "\n";
	std::cout << "      amounts to:        " << (sizeof fciD)*fciD.size()*fciD.size()*fciD.size()*fciD.size()/(1000.0*1000.0) << " MB\n";
	std::cout << "\n" << std::flush;

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


	//IFDBG std::cout << "\n";
	std::cout << "      first transformation..." << std::flush;
	// first transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<spinorSize; i++) {
		for (int nu=0; nu<spinorSize; nu++) {
			for (int lam=0; lam<spinorSize; lam++) {
				for (int sig=0; sig<spinorSize; sig++) {
					tmp[i][nu][lam][sig] = 0;
					for (int mu=0; mu<spinorSize; mu++) {
						tmp[i][nu][lam][sig] += spinor(mu, i) * fciD[mu][nu][lam][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;

	std::cout << "      second transformation..." << std::flush;
	// second transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			for (int lam=0; lam<spinorSize; lam++) {
				for (int sig=0; sig<spinorSize; sig++) {
					tmp2[i][j][lam][sig] = 0;
					for (int nu=0; nu<spinorSize; nu++) {
						tmp2[i][j][lam][sig] += std::conj(spinor(nu, j)) * tmp[i][nu][lam][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;

	std::cout << "      third transformation..." << std::flush;
	// third transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			for (int k=0; k<spinorSize; k++) {
				for (int sig=0; sig<spinorSize; sig++) {
					tmp3[i][j][k][sig] = 0;
					for (int lam=0; lam<spinorSize; lam++) {
						tmp3[i][j][k][sig] += spinor(lam, k) * tmp2[i][j][lam][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;

	std::cout << "      fourth transformation..." << std::flush;
	// fourth transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<spinorSize; i++) {
		for (int j=0; j<spinorSize; j++) {
			for (int k=0; k<spinorSize; k++) {
				for (int l=0; l<spinorSize; l++) {
					fourCenterIntegral[i][j][k][l] = 0;
					for (int sig=0; sig<spinorSize; sig++) {
						fourCenterIntegral[i][j][k][l] += std::conj(spinor(sig, l)) * tmp3[i][j][k][sig];
					}
				}
			}
		}
	}
	std::cout << " done.\n" << std::flush;
	std::cout << "\n" << std::flush;
	
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
	auto antisym = [&](int i, int j, int k, int l) { return fourCenterIntegral[i][k][j][l] - fourCenterIntegral[i][l][j][k]; };

	std::cout << "\n:: constructing A and B matricies...\n" << std::flush;

	A.resize( nocc*nvirt, nocc*nvirt );
	B.resize( nocc*nvirt, nocc*nvirt );
	std::cout << "      size of matricies:   2 x " << (nocc*nvirt) << " x " << (nocc*nvirt) << "\n";
	std::cout << "      amounts to:          " << 2.0*(sizeof A)*A.rows()*A.rows()/(1000.0*1000.0) << " MB\n";

	for (int i=0; i<nocc; i++) {
		for (int a=nocc; a<spinorSize; a++) {
			const int x = i*nvirt + a - nocc;

			for (int j=0; j<nocc; j++) {
				for (int b=nocc; b<spinorSize; b++) {
					//const int y = j*nvirt + b - nocc;
					//if (y > x) {
					//	continue;
					//}

					// (rb||as) = (rb|as) - (rb|sa)
					A( i*nvirt+a-nocc, j*nvirt+b-nocc ) = antisym(a, j, i, b);
					//A(x, y) = antisym(a, j, i, b);
					//A(y, x) = std::conj(antisym(a, j, i, b));
					if (a==b and i==j) {
						A( i*nvirt+a-nocc, j*nvirt+b-nocc ) += (epsilon[a] - epsilon[i]);
						//A(x, y) += (epsilon[a] - epsilon[i]);
					}
					// (rs||ab) = (rs|ab) - (rs|ba)
					B( i*nvirt+a-nocc, j*nvirt+b-nocc ) = antisym(a, b, i, j);
					//B(x, y) = antisym(a, b, i, j);
					//B(y, x) = antisym(a, b, i, j);
				}
			}
		}
	}

	// write A and B to files polly_a.re, polly_a.im, polly_b.re and polly_b.im
    std::cout << "      writing logfiles...\n" << std::flush;
	std::ofstream are, aim, bre, bim;
	are.open("polly_a.r");
	aim.open("polly_a.i");
	bre.open("polly_b.r");
	bim.open("polly_b.i");
	const int dim = A.rows()*(A.rows()-1)/2 + A.rows();
	are << std::setprecision(12) << dim << "\n";
	aim << std::setprecision(12) << dim << "\n";
	bre << std::setprecision(12) << dim << "\n";
	bim << std::setprecision(12) << dim << "\n";
	for (int x=0; x<A.rows(); x++) {
		for (int y=0; y<=x; y++) {
			{are << A(y, x).real() << "\n";}
			{aim << A(y, x).imag() << "\n";}
			{bre << B(y, x).real() << "\n";}
			{bim << B(y, x).imag() << "\n";}
		}
	}
	are.close();
	aim.close();
	bre.close();
	bim.close();

    // writing to control file
    std::ofstream cont;
    cont.open("control", std::ios_base::app);
    cont << "$polly_a_real        file=polly_a.r\n";
    cont << "$polly_a_imag        file=polly_a.i\n";
    cont << "$polly_b_real        file=polly_b.r\n";
    cont << "$polly_b_imag        file=polly_b.i\n";
    cont.close();

	std::cout << "\n" << std::flush;
}