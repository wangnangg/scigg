#pragma once
#include "matrix.hpp"
#include "spmatrix.hpp"
#include "vector.hpp"

namespace markovgg
{
// y = alpha * A * x + beta * y
void spblas_matrix_vector(real_t alpha, const spmatrix& A, bool transposeA,
                          vector_const_view x, real_t beta,
                          vector_mutable_view y);

// C = alpha * A * B + beta * C
void spblas_matrix_matrix(real_t alpha, const spmatrix& A, bool transposeA,
                          const spmatrix& B, bool transposeB, real_t beta,
                          spmatrix& C);
}
