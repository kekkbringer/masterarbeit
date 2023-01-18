#pragma once

#include <complex>
#include <vector>
#include <string>

typedef std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> fourD;

//void readFourcenter(std::complex<double> ****fci);
//template <int n>
//void readFourcenter(std::complex<double> (&fci)[n][n][n][n]);
//template <typename fourD>
//void readFourcenter(fourD &fci);
fourD readFourcenter(std::string location = "");
