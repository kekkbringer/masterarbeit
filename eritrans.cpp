#include "eritrans.hpp"
#include "read_fourcenter.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <string>
#include <complex>
#include <iomanip>

#include <Eigen/Core>

Eigen::VectorXcd& eritrans(const Eigen::MatrixXcd &spinorCAO, const int nocc, const int nvirt,
		Eigen::MatrixXcd& A, Eigen::MatrixXcd& B) {
	std::cout << " :: starting spinor transformation of fourcenter integrals\n" << std::flush;

	// reading coord file to know which atoms occur in which order
	//std::cout << "fourcenter reading coord...\n";
	std::ifstream coord("coord");
	if (coord.fail()) std::cout << "\nWARNING: could not read coord file!\n";
	std::string line;
	getline(coord, line); // $coord
	std::vector<std::string> atoms;
	while(getline(coord, line)) {
		if (line[0] == '#') continue;
		if (line[0] == '$') break;
		std::string word;
		std::istringstream iss(line);
		while(iss >> word) {} // now word contains the elemental symbol
		atoms.push_back(word);
	}
	coord.close();

	const int atomNum = atoms.size();

	// reading basis set info from 'basis' file
	std::ifstream basis("basis");
	if (basis.fail()) std::cout << "\nWARNING: could not read basis file!\n";

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
			std::istringstream bss(line);
			std::string word, a;
			bss >> word; // '#'
			bss >> a; // 'li'
			while (bss >> word) {
				if (word[0] == '[') {
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
							i++;
							g.push_back(0);
							for (; i<word.length()-1; i++) { // g
								if (word[i] == 'g') break;
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

	const int fciSize = sSize + 3*pSize + 6*dSize + 10*fSize + 15*gSize;	// nCAO
	const int returnSize = sSize + 3*pSize + 5*dSize + 7*fSize + 9*gSize;	// nSAO

	std::vector<std::vector<int>> orbitalIndexMap(sSize+pSize+dSize+fSize+gSize+1);
	orbitalIndexMap[0] = {-1};

	std::vector<std::vector<int>> returnMap(sSize+pSize+dSize+fSize+gSize+1);
	returnMap[0] = {-1};

	// first: all s orbitals
	int orbIndex = 1;
	int aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		//iterate over s-functions of current atom
		for (int i=0; i<s[atom]; i++) {
			orbitalIndexMap[orbIndex] = {aoIndex};
			orbIndex++;
			aoIndex++;
		}
		aoIndex += 3*p[atom] + 6*d[atom] + 10*f[atom] + 15*g[atom];
	}
	// second: all p orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]; // shift by s-functions of this atom
		// iterate over p-functions of current atom#
		for (int i=0; i<p[atom]; i++) {
			orbitalIndexMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2};
			orbIndex++;
			aoIndex += 3;
		}
		aoIndex += 6*d[atom] + 10*f[atom] + 15*g[atom];
	}
	// third: d orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]+3*p[atom]; // shift by s and p
		// iterate over d-functions of current atom
		for (int i=0; i<d[atom]; i++) {
			orbitalIndexMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2, aoIndex+3, aoIndex+4, aoIndex+5};
			orbIndex++;
			aoIndex += 6;
		}
		aoIndex += 10*f[atom] + 15*g[atom];
	}
	// fourth: f orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]+3*p[atom]+6*d[atom]; // shift by s, p and d
		// iterate over f-functions of current atom
		for (int i=0; i<f[atom]; i++) {
			orbitalIndexMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2, aoIndex+3, aoIndex+4,
						     aoIndex+5, aoIndex+6, aoIndex+7, aoIndex+8, aoIndex+9};
			orbIndex++;
			aoIndex += 10;
		}
		aoIndex += 15*g[atom];
	}
	// fifth: g-orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]+3*p[atom]+6*d[atom]+10*f[atom]; // shift by s, p, d and f
		// iterate over g-functions of current atom
		for (int i=0; i<g[atom]; i++) {
			orbitalIndexMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2, aoIndex+3, aoIndex+4,
						     aoIndex+5, aoIndex+6, aoIndex+7, aoIndex+8, aoIndex+9,
						     aoIndex+10, aoIndex+11, aoIndex+12, aoIndex+13, aoIndex+14};
			orbIndex++;
			aoIndex += 15;
		}
	}

	/**********************************************************************
	 *                            returnMap                               *
	 **********************************************************************
	 *
	 * TODO: comment
	 */

	// first: all s orbitals
	orbIndex = 1;
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		//iterate over s-functions of current atom
		for (int i=0; i<s[atom]; i++) {
			returnMap[orbIndex] = {aoIndex};
			orbIndex++;
			aoIndex++;
		}
		aoIndex += 3*p[atom] + 5*d[atom] + 7*f[atom] + 9*g[atom];
	}
	// second: all p orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]; // shift by s-functions of this atom
		// iterate over p-functions of current atom
		for (int i=0; i<p[atom]; i++) {
			returnMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2};
			orbIndex++;
			aoIndex += 3;
		}
		aoIndex += 5*d[atom] + 7*f[atom] + 9*g[atom];
	}
	// third: d orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]+3*p[atom]; // shift by s and p
		// iterate over d-functions of current atom
		for (int i=0; i<d[atom]; i++) {
			returnMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2, aoIndex+3, aoIndex+4};
			orbIndex++;
			aoIndex += 5;
		}
		aoIndex += 7*f[atom] + 9*g[atom];
	}
	// fourth: f orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]+3*p[atom]+5*d[atom]; // shift by s, p and d
		// iterate over f-functions of current atom
		for (int i=0; i<f[atom]; i++) {
			returnMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2, aoIndex+3, aoIndex+4, aoIndex+5, aoIndex+6};
			orbIndex++;
			aoIndex += 7;
		}
		aoIndex += 9*g[atom];
	}
	// fifth: g-orbitals
	aoIndex = 0;
	for (int atom=0; atom<atomNum; atom++) {
		aoIndex += s[atom]+3*p[atom]+5*d[atom]+7*f[atom]; // shift by s, p, d and f
		// iterate over g-functions of current atom
		for (int i=0; i<g[atom]; i++) {
			returnMap[orbIndex] = {aoIndex, aoIndex+1, aoIndex+2, aoIndex+3, aoIndex+4, aoIndex+5, aoIndex+6, aoIndex+7, aoIndex+8};
			orbIndex++;
			aoIndex += 9;
		}
	}

	Eigen::MatrixXcd swapMat = Eigen::MatrixXcd::Zero(sSize+3*pSize+6*dSize+10*fSize+15*gSize, sSize+3*pSize+6*dSize+10*fSize+15*gSize);
	int orbcounter = 0;
	for (int i=1; i<sSize+pSize+dSize+fSize+gSize+1; i++) {
		for (auto x: orbitalIndexMap[i]) {
			swapMat(orbcounter, x) = 1;
			orbcounter++;
		}
	}

	// save swap matrix to file berryswap.r
	std::ofstream swapfile;
	swapfile.open("berryswap.r");
	swapfile << sSize+3*pSize+6*dSize+10*fSize+15*gSize << "\n";
	for (int i=0; i<sSize+3*pSize+6*dSize+10*fSize+15*gSize; i++) {
		for (int j=0; j<sSize+3*pSize+6*dSize+10*fSize+15*gSize; j++) {
			swapfile << std::fixed << std::setprecision(15) << swapMat(i, j).real() << "\n";
		}
	}
	swapfile.close();


	int j=0;



	// init A and B matrices
	A.resize( nocc*nvirt, nocc*nvirt );
	B.resize( nocc*nvirt, nocc*nvirt );

	static Eigen::VectorXcd ailkasym(nvirt*nocc*nocc*nocc);


	// init Four Center Integrals
	const int nCAO = fciSize;
	const int nSAO = returnSize;
	std::cout << "eritrafo memory requirement for first trafo step: " << sizeof(std::complex<double>)*2.0*nCAO*nCAO*nCAO/1000.0/1000.0 << " MB\n";
	int batchnum;

	for (batchnum=nocc; batchnum<nocc+nvirt; batchnum++) {
		// ============================================== beginning of one batch =======================================
		// reading acutal real part
		std::ifstream re("fourcenter.r");
		std::ifstream im("fourcenter.i");
		std::string lineIm;

		std::cout << "   " << batchnum-nocc+1 << "/" << nvirt << std::endl;

		std::vector<Eigen::MatrixXcd> trans1(2*nCAO, Eigen::MatrixXcd::Zero(nCAO, nCAO));

		int i, a, b; // [ij||ab]
		while (getline(re, line)) {
			getline(im, lineIm);
			// determine indicies and degeneracy
			std::istringstream iss(line);
			std::string word;
			int degeneracy = 1;
			iss >> word; i = stoi(word);
			iss >> word; j = stoi(word);
			iss >> word; a = stoi(word);
			getline(re, line);
			getline(im, lineIm);
			b = stoi(line);

			// mega unclean aber scheiß drauf
			const unsigned int id1 = i + 1'0000*j + 1'0000'0000*a + 1'0000'0000'0000*b;
			const unsigned int id2 = a + 1'0000*b + 1'0000'0000*i + 1'0000'0000'0000*j;
			const unsigned int id3 = j + 1'0000*i + 1'0000'0000*b + 1'0000'0000'0000*a;
			const unsigned int id4 = b + 1'0000*a + 1'0000'0000*j + 1'0000'0000'0000*i;

			// reading in integerals
			for (const auto& alpha: orbitalIndexMap.at(i)) {
				for (const auto& beta: orbitalIndexMap.at(j)) {
					for (const auto& gamma: orbitalIndexMap.at(a)) {
						for (const auto& delta: orbitalIndexMap.at(b)) {
							getline(re, line);
							getline(im, lineIm);
							const double valRe = stod(line);
							const double valIm = stod(lineIm);

							// spin up + spin up
							trans1[beta](gamma, delta)
								+= spinorCAO(alpha, batchnum)*std::complex<double>(valRe,  valIm);
							if (id1!=id2)
							trans1[delta](alpha, beta)
								+= spinorCAO(gamma, batchnum)*std::complex<double>(valRe,  valIm);
							if (id1!=id3 and id2!=id3)
							trans1[alpha](delta, gamma)
								+= spinorCAO( beta, batchnum)*std::complex<double>(valRe, -valIm);
							if (id1!=id4 and id2!=id4 and id3!=id4)
							trans1[gamma](beta, alpha)
								+= spinorCAO(delta, batchnum)*std::complex<double>(valRe, -valIm);

							// spin down + spin down
							trans1[beta+nCAO](gamma, delta)
								+= spinorCAO(alpha+nCAO, batchnum)*std::complex<double>(valRe,  valIm);
							if (id1!=id2)
							trans1[delta+nCAO](alpha, beta)
								+= spinorCAO(gamma+nCAO, batchnum)*std::complex<double>(valRe,  valIm);
							if (id1!=id3 and id2!=id3)
							trans1[alpha+nCAO](delta, gamma)
								+= spinorCAO( beta+nCAO, batchnum)*std::complex<double>(valRe, -valIm);
							if (id1!=id4 and id2!=id4 and id3!=id4)
							trans1[gamma+nCAO](beta, alpha)
								+= spinorCAO(delta+nCAO, batchnum)*std::complex<double>(valRe, -valIm);
						}
					}
				}
			}
		}
		re.close();
		im.close();


		// second transformation step
		// (a nu | kap lam)  ->  (a p | kap lam)
		std::vector<Eigen::MatrixXcd> trans2(2*nCAO, Eigen::MatrixXcd::Zero(nCAO, nCAO));
		for (int p=0; p<2*nSAO; p++) {
			for (int kap=0; kap<nCAO; kap++) {
				for (int lam=0; lam<nCAO; lam++) {
					for (int nu=0; nu<2*nCAO; nu++) {
						trans2[p](kap, lam) += std::conj(spinorCAO(nu, p)) * trans1[nu](kap, lam);
					}
				}
			}
		}
		trans1.resize(0);
		trans1.shrink_to_fit();
		
		// third transformation step
		// (a p | kap lam)  ->  (a p | q lam)
		std::vector<Eigen::MatrixXcd> trans3(2*nCAO, Eigen::MatrixXcd::Zero(2*nCAO, 2*nCAO));
		for (int p=0; p<2*nSAO; p++) {
			for (int q=0; q<2*nSAO; q++) {
				for (int lam=0; lam<nCAO; lam++) {
					for (int kap=0; kap<nCAO; kap++) {
						// spin up + spin up
						trans3[p](q, lam) += spinorCAO(kap, q) * trans2[p](kap, lam);
						// spin down + spin down
						trans3[p](q, lam+nCAO) += spinorCAO(kap+nCAO, q) * trans2[p](kap, lam);
					}
				}
			}
		}
		trans2.resize(0);
		trans2.shrink_to_fit();


		// fourth transformation
		// (a p | q lam)  ->  (a p | q r)
		std::vector<Eigen::MatrixXcd> trans4(2*nCAO, Eigen::MatrixXcd::Zero(2*nCAO, 2*nCAO));
		for (int p=0; p<2*nSAO; p++) {
			for (int q=0; q<2*nSAO; q++) {
				for (int r=0; r<2*nSAO; r++) {
					for (int lam=0; lam<2*nCAO; lam++) {
						trans4[p](q, r) += std::conj(spinorCAO(lam, r)) * trans3[p](q, lam);
					}
				}
			}
		}
		trans3.resize(0);
		trans3.shrink_to_fit();

		
		// write to A and B
		for (int i=0; i<nocc; i++) {
			for (int j=0; j<nocc; j++) {
				for (int b=nocc; b<nocc+nvirt; b++) {
					const int a = batchnum;
					A ( i*nvirt+a-nocc, j*nvirt+b-nocc ) = trans4[i](j, b) - trans4[b](j, i);
					B ( i*nvirt+a-nocc, j*nvirt+b-nocc ) = trans4[i](b, j) - trans4[j](b, i);
				}
				
				for (int s=0; s<nocc; s++) { // a i j s
					const int a = batchnum;
					const int index = (a-nocc) + nvirt*s + nvirt*nocc*j + nvirt*nocc*nocc*i;
					const auto tmp = std::conj(trans4[i](s, j) - trans4[j](s, i));
					ailkasym(index) = tmp;
				}
			}
		}
		// ============================================== end of one batch =============================================
	}
	std::cout << "eritrans ended successfully" << std::endl;

	return ailkasym;
}
