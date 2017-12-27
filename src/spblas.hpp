#pragma once
#include "matrix.hpp"
#include "spmatrix.hpp"
#include "vector.hpp"

namespace markovgg
{
// y = alpha * A * x + beta * y
void spblas_matrix_vector(real_t alpha, spmatrix_const_view A,
                          vector_const_view x, real_t beta,
                          vector_mutable_view y);

// C = alpha * A * B + beta * C
}
