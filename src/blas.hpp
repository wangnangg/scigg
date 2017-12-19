#pragma once

#include <cassert>
#include "matrix.hpp"
#include "vector.hpp"
extern "C" {
#include "cblas.h"
}
namespace markovgg
{
// level 1 blas

real_t blas_dot(vector_const_view v1, vector_const_view v2);

real_t blas_norm2(vector_const_view v);

real_t blas_abs_sum(vector_const_view v);

size_t blas_abs_max_idx(vector_const_view v);

void blas_swap(vector_mutable_view v1, vector_mutable_view v2);

void blas_copy(vector_const_view src, vector_mutable_view dst);

// y = ax + y
void blas_axpy(real_t a, vector_const_view x, vector_mutable_view y);

// Rotation, omitted
// void cblas_drotg(double *a, double *b, double *c, double *s);
// void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double
// *P);  void cblas_drot(const int N, double *X, const int incX,
//               double *Y, const int incY, const double c, const double  s);
// void cblas_drotm(const int N, double *X, const int incX,
//                double *Y, const int incY, const double *P);

// x = ax
void blas_scale(real_t alpha, vector_mutable_view v);

// level 2 blas
// y = alpha * A * x + beta * y
void blas_matrix_vector(real_t alpha, matrix_const_view A, bool transposeA,
                        vector_const_view x, real_t beta,
                        vector_mutable_view y);

// rank 1 or 2 op, omitted
// A = alpha * x * y' + A
// A = alpha * x * y' + beta * y * x' + A

// level 3 blas
// C = alpha * A * B + beta * C
void blas_matrix_matrix(real_t alpha, matrix_const_view A, bool transposeA,
                        matrix_const_view B, bool transposeB, real_t beta,
                        matrix_mutable_view C);

// rank k or 2k op, omitted
};
