#pragma once

#include <Eigen/Core>

#include "read_fourcenter.hpp"

Eigen::VectorXcd berryRHS(const int nuc, const int cart, const Eigen::VectorXcd& ailkasym);
