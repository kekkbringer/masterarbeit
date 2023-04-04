#pragma once

#include <Eigen/Core>

#include "read_fourcenter.hpp"

void split1efiles(int atomnum);
Eigen::VectorXcd berryRHS(const int nuc, const int cart, const Eigen::VectorXcd& ailkasym);
