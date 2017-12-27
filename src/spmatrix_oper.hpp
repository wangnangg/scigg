#pragma once

#include "spmatrix.hpp"
#include "vector.hpp"

namespace markovgg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         bool is_row_compressed);

// y = A . x
void dot(vector_mutable_view y, spmatrix_const_view A, vector_const_view x);

vector dot(spmatrix_const_view A, vector_const_view x);
}
