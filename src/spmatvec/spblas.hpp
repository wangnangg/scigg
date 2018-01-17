#pragma once
#include "matvec.hpp"
#include "spmatrix.hpp"

namespace scigg
{
// y = ax + y
void spblas_zero_fillin_axpy(real_t alpha, spvec_const_view x,
                             spvec_mutable_view y);
void spblas_axpy(real_t alpha, spvec_const_view x, vector_mutable_view y);

// y = alpha * A * x + beta * y
void spblas_matrix_vector(real_t alpha, spmatrix_const_view A,
                          vector_const_view x, real_t beta,
                          vector_mutable_view y);

// C = A + B
void spblas_matrix_add(spmatrix& C, spmatrix_const_view A,
                       spmatrix_const_view B);

// C = alpha * A * B
void spblas_matrix_matrix(spmatrix& C, real_t alpha, spmatrix_const_view A,
                          spmatrix_const_view B);
}
