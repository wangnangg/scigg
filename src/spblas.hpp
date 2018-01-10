#pragma once
#include "matrix.hpp"
#include "spmatrix.hpp"
#include "vector.hpp"

namespace markovgg
{
real_t dot(spvec_const_view x, vector_const_view y);
real_t dot(spvec_const_view x, spvec_const_view y);
inline real_t dot(vector_const_view x, spvec_const_view y) { return dot(y, x); }
inline real_t operator*(spvec_const_view x, vector_const_view y)
{
    return dot(x, y);
}
inline real_t operator*(vector_const_view x, spvec_const_view y)
{
    return dot(y, x);
}
inline real_t operator*(spvec_const_view x, spvec_const_view y)
{
    return dot(x, y);
}

// y = ax + y
void spblas_zero_fillin_axpy(real_t alpha, spvec_const_view x,
                             spvec_mutable_view y);
void spblas_axpy(real_t alpha, spvec_const_view x, vector_mutable_view y);

// y = alpha * A * x + beta * y
void spblas_matrix_vector(real_t alpha, spmatrix_const_view A,
                          vector_const_view x, real_t beta,
                          vector_mutable_view y);

// C = alpha * A * B + beta * C
}
