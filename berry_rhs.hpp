#pragma once

#include <Eigen/Core>

#include "read_fourcenter.hpp"

void splitExchange(int atomNum, int ncao);
void splitJSxi(int atomNum, int ncao);
void split1efiles(int atomnum, int ncao);
Eigen::VectorXcd berryRHS(const int nuc, const int cart, const Eigen::VectorXcd& ailkasym);
