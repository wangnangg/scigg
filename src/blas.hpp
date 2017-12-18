#pragma once

#include "matrix.hpp"
#include "vector.hpp"

namespace markovgg
{
// level 1 blas

real_t dot(const vector& v1, const vector& v2);

real_t norm2(const vector& v);

real_t abs_sum(const vector& v);

size_t abs_max_idx(const vector& v);

void swap(vector& v1, vector& v2);

void copy(const vector& src, vector& dst);

// y = ax + y
void axpy(real_t a, const vector& x, vector& y);

// Rotation, omitted
// void cblas_drotg(double *a, double *b, double *c, double *s);
// void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double
// *P);  void cblas_drot(const int N, double *X, const int incX,
//               double *Y, const int incY, const double c, const double  s);
// void cblas_drotm(const int N, double *X, const int incX,
//                double *Y, const int incY, const double *P);

// x = ax
void scale(real_t alpha, vector& v);

// level 2 blas

// y = alpha * A * x + beta * y
void matrix_vector(real_t alpha, const matrix& A, bool transpose,
                   const vector& x, real_t beta, vector& y);

// rank 1 or 2 op, omitted
// A = alpha * x * y' + A
// A = alpha * x * y' + beta * y * x' + A

// level 3 blas
// C = alpha * A * B + beta * C
void gen_matrix_matrix(real_t alpha, const matrix& A, bool tranposeA,
                       const matrix& B, bool transposeB, real_t beta,
                       matrix& C);

// rank k or 2k op, omitted
};
