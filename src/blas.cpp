#include "blas.hpp"
#include <utility>

namespace markovgg
{
// level 1 blas

real_t blas_dot(vector_const_view v1, vector_const_view v2)
{
    assert(v1.dim() == v2.dim());
    return cblas_ddot(v1.dim(), &v1[0], v1.inc(), &v2[0], v2.inc());
}
real_t blas_norm2(vector_const_view v)
{
    return cblas_dnrm2(v.dim(), &v[0], v.inc());
}

real_t blas_abs_sum(vector_const_view v)
{
    return cblas_dasum(v.dim(), &v[0], v.inc());
}

size_t blas_abs_max_idx(vector_const_view v)
{
    return cblas_idamax(v.dim(), &v[0], v.inc());
}

void blas_swap(vector_mutable_view v1, vector_mutable_view v2)
{
    assert(v1.dim() == v2.dim());
    cblas_dswap(v1.dim(), &v1[0], v1.inc(), &v2[0], v2.inc());
}

void blas_copy(vector_const_view src, vector_mutable_view dst)
{
    assert(src.dim() == dst.dim());
    cblas_dcopy(src.dim(), &src[0], src.inc(), &dst[0], dst.inc());
}

// y = ax + y
void blas_axpy(real_t a, vector_const_view x, vector_mutable_view y)
{
    assert(x.dim() == y.dim());
    cblas_daxpy(x.dim(), a, &x[0], x.inc(), &y[0], y.inc());
}

// Rotation, omitted
// void cblas_drotg(double *a, double *b, double *c, double *s);
// void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double
// *P);  void cblas_drot(const int N, double *X, const int incX,
//               double *Y, const int incY, const double c, const double  s);
// void cblas_drotm(const int N, double *X, const int incX,
//                double *Y, const int incY, const double *P);

// x = ax
void blas_scale(real_t alpha, vector_mutable_view v)
{
    cblas_dscal(v.dim(), alpha, &v[0], v.inc());
}

// level 2 blas

// y = alpha * A * x + beta * y
void blas_matrix_vector(real_t alpha, matrix_const_view A, bool transposeA,
                        vector_const_view x, real_t beta, vector_mutable_view y)
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
                A.n(), alpha, &A(0, 0), A.ldim(), &x[0], x.inc(), beta, &y[0],
                y.inc());
}

// A = alpha*x*y' +  A
void blas_rank1(real_t alpha, vector_const_view x, vector_const_view y,
                matrix_mutable_view A)
{
    cblas_dger(CblasRowMajor, A.m(), A.n(), alpha, &x[0], x.inc(), &y[0],
               y.inc(), &A(0, 0), A.ldim());
}

// rank 1 or 2 op, omitted
// A = alpha * x * y' + A
// A = alpha * x * y' + beta * y * x' + A

// level 3 blas
// C = alpha * A * B + beta * C
void blas_matrix_matrix(real_t alpha, matrix_const_view A, bool transposeA,
                        matrix_const_view B, bool transposeB, real_t beta,
                        matrix_mutable_view C)
{
    size_t K = transposeA ? A.m() : A.n();
    assert(transposeB ? B.n() : B.m() == K);
    assert(transposeA ? A.n() : A.m() == C.m());
    assert(transposeB ? B.m() : B.n() == C.n());

    cblas_dgemm(CblasRowMajor, transposeA ? CblasTrans : CblasNoTrans,
                transposeB ? CblasTrans : CblasNoTrans, C.m(), C.n(), K, alpha,
                &A(0, 0), A.ldim(), &B(0, 0), B.ldim(), beta, &C(0, 0),
                C.ldim());
}

// rank k or 2k op, omitted
};
