#pragma once
#include <utility>
#include <vector>
#include "type.hpp"
#include "matrix.hpp"

using namespace markovgg;

void print_sep();
vec create_vec(const std::vector<real_t>& val);
sqr_mat create_sqr_mat(size_t dim, const std::vector<real_t>& val);
matrix create_matrix(size_t row_dim, size_t col_dim,
                     const std::vector<real_t>& val);
