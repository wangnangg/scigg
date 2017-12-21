#pragma once

#include "spmatrix.hpp"
#include "vector.hpp"

namespace markovgg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         spmatrix_format format);

// y = A . x
void dot(vector_mutable_view y, const spmatrix& A, bool tranposeA,
         vector_const_view x);

vector dot(const spmatrix& A, bool transposeA, vector_const_view x);
}
