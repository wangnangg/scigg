#include "blas.hpp"
#include <cassert>
#include <utility>
extern "C" {
#include "cblas.h"
}
namespace markovgg
{
// level 1 blas

real_t dot(const vector& v1, const vector& v2)
{
    return cblas_ddot(v1.dim(), v1.begin(), 1, v2.begin(), 1);
}

real_t norm2(const vector& v) { return cblas_dnrm2(v.dim(), v.begin(), 1); }

real_t abs_sum(const vector& v) { return cblas_dasum(v.dim(), v.begin(), 1); }

size_t abs_max_idx(const vector& v)
{
    return cblas_idamax(v.dim(), v.begin(), 1);
}

void swap(vector& v1, vector& v2)
{
    assert(v1.dim() == v2.dim());
    std::swap(v1, v2);
}

void copy(const vector& src, vector& dst)
{
    assert(src.dim() == dst.dim());
    cblas_dcopy(src.dim(), src.begin(), 1, dst.begin(), 1);
}

// y = ax + y
void axpy(real_t a, const vector& x, vector& y)
{
    assert(x.dim() == y.dim());
    cblas_daxpy(x.dim(), a, x.begin(), 1, y.begin(), 1);
}

// Rotation, omitted
// void cblas_drotg(double *a, double *b, double *c, double *s);
// void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double
// *P);  void cblas_drot(const int N, double *X, const int incX,
//               double *Y, const int incY, const double c, const double  s);
// void cblas_drotm(const int N, double *X, const int incX,
//                double *Y, const int incY, const double *P);

// x = ax
void scale(real_t alpha, vector& v)
{
    cblas_dscal(v.dim(), alpha, v.begin(), 1);
}

// level 2 blas

// y = alpha * A * x + beta * y
void matrix_vector(real_t alpha, const matrix& A, bool transposeA,
                   const vector& x, real_t beta, vector& y)
{
    if (!transposeA)
    {
        assert(A.n() == x.dim());
        assert(A.m() == y.dim());
    }
    else
    {
        assert(A.m() == x.dim());
        assert(A.n() == y.dim());
    }
    cblas_dgemv(CblasRowMajor, transposeA ? CblasTrans : CblasNoTrans, A.m(),
                A.n(), alpha, A.begin(), A.n(), x.begin(), 1, beta, y.begin(),
                1);
}

// rank 1 or 2 op, omitted
// A = alpha * x * y' + A
// A = alpha * x * y' + beta * y * x' + A

// level 3 blas
// C = alpha * A * B + beta * C
void gen_matrix_matrix(real_t alpha, const matrix& A, bool transposeA,
                       const matrix& B, bool transposeB, real_t beta, matrix& C)
{
    size_t K = transposeA ? A.m() : A.n();
    assert(transposeB ? B.n() : B.m() == K);
    assert(transposeA ? A.n() : A.m() == C.m());
    assert(transposeB ? B.m() : B.n() == C.n());

    cblas_dgemm(CblasRowMajor, transposeA ? CblasTrans : CblasNoTrans,
                transposeB ? CblasTrans : CblasNoTrans, C.m(), C.n(), K, alpha,
                A.begin(), A.n(), B.begin(), B.n(), beta, C.begin(), C.n());
}

// rank k or 2k op, omitted
};
