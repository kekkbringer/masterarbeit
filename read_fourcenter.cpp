#include "read_fourcenter.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <complex>
#include <iomanip>

fourD readFourcenter() {
	//std::cout << "\n\n========================================\nreading four center integrals... oh boi...\n";

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
		//std::cout << "elemental symbol: " << word << "\n";
		atoms.push_back(word);
	}
	coord.close();

	const int atomNum = atoms.size();
	//std::cout << "\nnumber of atoms: " << atomNum << "\n";



	// reading basis set info from 'basis' file
	std::ifstream basis("basis");
	if (basis.fail()) std::cout << "\nWARNING: could not read basis file!\n";
	std::cout << "      minor warning: I'm reading from comments in basis file, e.g. '# li  (7s3p) / [3s2p]    {511/21}'\n";
	std::cout << "      minor warning: only s-, p- and d-type orbitals are supported yet!\n";

	std::vector<int> s(atomNum);
	std::vector<int> p(atomNum);
	std::vector<int> d(atomNum);
	int sSize = 0;
	int pSize = 0;
	int dSize = 0;
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
						}
					}
				}
			}
		}
	}

	const int fciSize = sSize + 3*pSize + 6*dSize;
	const int returnSize = sSize + 3*pSize + 5*dSize;

	std::vector<std::vector<int>> orbitalIndexMap(sSize+pSize+dSize+1);
	orbitalIndexMap[0] = {-1};

	std::vector<std::vector<int>> returnMap(sSize+pSize+dSize+1);
	returnMap[0] = {-1};

	//std::cout << "init orbitalIndexMap...\n";
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
		aoIndex += 3*p[atom] + 6*d[atom];
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
		aoIndex += 6*d[atom];
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
	}

	/**********************************************************************
	 *                            returnMap                               *
	 **********************************************************************
	 *
	 * TODO: comment
	 */

	//std::cout << "init return Map...\n";
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
		aoIndex += 3*p[atom] + 5*d[atom];
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
		aoIndex += 5*d[atom];
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
	}
	int j=0;



	// init Four Center Integrals
	//std::cout << "size of four center integral: Tensor of rank 4 with dimension " << fciSize << "\n";
	fourD fci(fciSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(fciSize,
				std::vector<std::vector<std::complex<double>>>(fciSize,
					std::vector<std::complex<double>>(fciSize))));

	bool index[fciSize][fciSize][fciSize][fciSize];

	for (int i=0; i<fciSize; i++) {
		for (int j=0; j<fciSize; j++) {
			for (int k=0; k<fciSize; k++) {
				for (int l=0; l<fciSize; l++) {
					index[i][j][k][l] = false;
				}		
			}
		}
	}



	// reading acutal real part
	std::ifstream re("fourcenter.r");
	std::ifstream im("fourcenter.i");
	std::string lineIm;

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
		degeneracy *= orbitalIndexMap.at(i).size()
		            * orbitalIndexMap.at(j).size()
		            * orbitalIndexMap.at(a).size()
		            * orbitalIndexMap.at(b).size();
		//std::cout << "indicies: " << i << ", " << j << ", " << a << ", " << b << "\n";
		//std::cout << "	degeneracy: " << degeneracy << "\n";

		// reading in integerals
		for (const auto& alpha: orbitalIndexMap.at(i)) {
			for (const auto& beta: orbitalIndexMap.at(j)) {
				for (const auto& gamma: orbitalIndexMap.at(a)) {
					for (const auto& delta: orbitalIndexMap.at(b)) {
						getline(re, line);
						getline(im, lineIm);
						const double valRe = stod(line);
						const double valIm = stod(lineIm);

						fci[alpha][beta][gamma][delta] = std::complex<double>(valRe,  valIm);
						fci[gamma][delta][alpha][beta] = std::complex<double>(valRe,  valIm);
						fci[beta][alpha][delta][gamma] = std::complex<double>(valRe, -valIm);
						fci[delta][gamma][beta][alpha] = std::complex<double>(valRe, -valIm);

						index[alpha][beta][gamma][delta] = true;
						index[gamma][delta][alpha][beta] = true;
						index[beta][alpha][delta][gamma] = true;
						index[delta][gamma][beta][alpha] = true;
					}
				}
			}
		}
	}
	re.close();
	im.close();


	bool full = true;
	for (int i=0; i<fciSize; i++) {
		for (int j=0; j<fciSize; j++) {
			for (int k=0; k<fciSize; k++) {
				for (int l=1; l<fciSize; l++) {
					if(index[i][j][k][l] == false) {
						full = false;
						//std::cout << "indecies: " << i << ", " << j << ", " << k << ", " << l << "\n";
					}
				}		
			}
		}
	}
	if (full == false) {
		std::cout << "\n\nWARNING: fourcenterintegrals was not fully defined!\n\n";
	}

	// reading acutal imag part
	//std::ifstream im("fourcenter.i");

	//while (getline(im, line)) {
	//	// determine indicies and degeneracy
	//	std::istringstream iss(line);
	//	std::string word;
	//	int degeneracy = 1;
	//	iss >> word; i = stoi(word);
	//	iss >> word; j = stoi(word);
	//	iss >> word; a = stoi(word);
	//	getline(im, line);
	//	b = stoi(line);
	//	degeneracy *= orbitalIndexMap.at(i).size()
	//	            * orbitalIndexMap.at(j).size()
	//	            * orbitalIndexMap.at(a).size()
	//	            * orbitalIndexMap.at(b).size();
	//	//std::cout << "indicies: " << i << ", " << j << ", " << a << ", " << b << "\n";
	//	//std::cout << "	degeneracy: " << degeneracy << "\n";

	//	// reading in integerals
	//	for (const auto& alpha: orbitalIndexMap.at(i)) {
	//		for (const auto& beta: orbitalIndexMap.at(j)) {
	//			for (const auto& gamma: orbitalIndexMap.at(a)) {
	//				for (const auto& delta: orbitalIndexMap.at(b)) {
	//					getline(im, line);
	//					const double val = stod(line);
	//					fci[alpha][beta][gamma][delta] += std::complex<double>(0,  val);
	//					fci[gamma][delta][alpha][beta] += std::complex<double>(0,  val);
	//					fci[beta][alpha][delta][gamma] += std::complex<double>(0, -val);
	//					fci[delta][gamma][beta][alpha] += std::complex<double>(0, -val);
	//				}
	//			}
	//		}
	//	}
	//}
	//im.close();

	std::cout << "\n";
	//std::cout << "\n\ndone reading fcis\n=====================================================\n\n";
	


	// if no d orbitals are present we are done
	if (dSize == 0) {
		//std::cout << " no d-orbitals found, no transformation needed.\n";
		return fci;
	} // else we need to transform...
	


	/*****************************************************************************
	*         transform d orbitals from cartesian to real spherical              *
	*****************************************************************************/
	std::cout << " transforming d orbitals from cartesian to real spherical...\n";
	constexpr double normZZ   = 1.0/sqrt(12.0);//0.2886751346;//1.0/sqrt(6.0);
	constexpr double normXXYY = 0.5;//1.0/sqrt(2.0);
	constexpr double beta = 1.0;//0.655722;
	// calculate transformation matrix
	std::vector<std::vector<double> > trans(
			fciSize,
			std::vector<double>(returnSize));
	for (int orb=1; orb<=returnMap.size(); orb++) {
		// check degeneracy
		if (returnMap[orb].size() == 1) { // s
			trans[orbitalIndexMap[orb][0]][returnMap[orb][0]] = 1.0;
			//std::cout << orbitalIndexMap[orb][0] << " -> " << returnMap[orb][0] << "\n";
		} else if (returnMap[orb].size() == 3) { // p
			trans[orbitalIndexMap[orb][0]][returnMap[orb][0]] = 1.0;
			trans[orbitalIndexMap[orb][1]][returnMap[orb][1]] = 1.0;
			trans[orbitalIndexMap[orb][2]][returnMap[orb][2]] = 1.0;
			//std::cout << orbitalIndexMap[orb][0] << " -> " << returnMap[orb][0] << "\n";
			//std::cout << orbitalIndexMap[orb][1] << " -> " << returnMap[orb][1] << "\n";
			//std::cout << orbitalIndexMap[orb][2] << " -> " << returnMap[orb][2] << "\n";
		} else if (returnMap[orb].size() == 5) { // d
			// dz^2
			trans[orbitalIndexMap[orb][0]][returnMap[orb][0]] = -1.0 * normZZ * beta;
			trans[orbitalIndexMap[orb][1]][returnMap[orb][0]] = -1.0 * normZZ * beta;
			trans[orbitalIndexMap[orb][2]][returnMap[orb][0]] =  2.0 * normZZ * beta;
			//std::cout << "2*" << orbitalIndexMap[orb][2] << " - " << orbitalIndexMap[orb][0] << " - " << orbitalIndexMap[orb][1] << " -> " << returnMap[orb][0] << "\n";
			// dxz
			trans[orbitalIndexMap[orb][4]][returnMap[orb][1]] = 1.0 * beta;
			//std::cout << orbitalIndexMap[orb][4] << " -> " << returnMap[orb][1] << "\n";
			// dyz
			trans[orbitalIndexMap[orb][5]][returnMap[orb][2]] = 1.0 * beta;
			//std::cout << orbitalIndexMap[orb][5] << " -> " << returnMap[orb][2] << "\n";
			// dxy
			trans[orbitalIndexMap[orb][3]][returnMap[orb][3]] = 1.0 * beta;
			//std::cout << orbitalIndexMap[orb][3] << " -> " << returnMap[orb][3] << "\n";
			// dx^2-y^2
			trans[orbitalIndexMap[orb][0]][returnMap[orb][4]] =  1.0 * normXXYY * beta;
			trans[orbitalIndexMap[orb][1]][returnMap[orb][4]] = -1.0 * normXXYY * beta;
			//std::cout << orbitalIndexMap[orb][0] << " - " << orbitalIndexMap[orb][1] << " -> " << returnMap[orb][4] << "\n";

		}
	}

	// print trans matrix
	//std::cout << std::setprecision(1);
	//for (auto x: trans) {
	//	for (auto y: x) {
	//		std::cout << y << "  ";
	//	}
	//	std::cout << "\n";
	//}
	//std::cout << "\n" << std::flush;


	// init return vector
	fourD returnfci(returnSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(returnSize,
				std::vector<std::vector<std::complex<double>>>(returnSize,
					std::vector<std::complex<double>>(returnSize))));

	fourD tmp(fciSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(fciSize,
				std::vector<std::vector<std::complex<double>>>(fciSize,
					std::vector<std::complex<double>>(fciSize))));
	fourD tmp2(fciSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(fciSize,
				std::vector<std::vector<std::complex<double>>>(fciSize,
					std::vector<std::complex<double>>(fciSize))));
	fourD tmp3(fciSize,
			std::vector<std::vector<std::vector<std::complex<double>>>>(fciSize,
				std::vector<std::vector<std::complex<double>>>(fciSize,
					std::vector<std::complex<double>>(fciSize))));


	// first transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<fciSize; i++) {
		for (int nu=0; nu<fciSize; nu++) {
			for (int lam=0; lam<fciSize; lam++) {
				for (int sig=0; sig<fciSize; sig++) {
					tmp[i][nu][lam][sig] = 0;
					for (int mu=0; mu<fciSize; mu++) {
						tmp[i][nu][lam][sig] += trans[mu][i] * fci[mu][nu][lam][sig];
					}
				}
			}
		}
	}

	// second transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<fciSize; i++) {
		for (int j=0; j<fciSize; j++) {
			for (int lam=0; lam<fciSize; lam++) {
				for (int sig=0; sig<fciSize; sig++) {
					tmp2[i][j][lam][sig] = 0;
					for (int nu=0; nu<fciSize; nu++) {
						tmp2[i][j][lam][sig] += trans[nu][j] * tmp[i][nu][lam][sig];
					}
				}
			}
		}
	}

	// third transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<fciSize; i++) {
		for (int j=0; j<fciSize; j++) {
			for (int k=0; k<fciSize; k++) {
				for (int sig=0; sig<fciSize; sig++) {
					tmp3[i][j][k][sig] = 0;
					for (int lam=0; lam<fciSize; lam++) {
						tmp3[i][j][k][sig] += trans[lam][k] * tmp2[i][j][lam][sig];
					}
				}
			}
		}
	}

	// fourth transformation
	#pragma omp parallel for collapse(4)
	for (int i=0; i<returnSize; i++) {
		for (int j=0; j<returnSize; j++) {
			for (int k=0; k<returnSize; k++) {
				for (int l=0; l<returnSize; l++) {
					returnfci[i][j][k][l] = 0;
					for (int sig=0; sig<fciSize; sig++) {
						returnfci[i][j][k][l] += trans[sig][l] * tmp3[i][j][k][sig];
					}
				}
			}
		}
	}



	return returnfci;
}
