#include "spblas.hpp"
#include "blas.hpp"
namespace markovgg
{
// y = alpha * A * x + beta * y
void spblas_matrix_vector(real_t alpha, const spmatrix& A, bool transposeA,
                          vector_const_view x, real_t beta,
                          vector_mutable_view y)
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
    blas_scale(beta, y);
    if ((!transposeA && A.format() == CPR_ROW) ||
        (transposeA && A.format() == CPR_COL))
    {
        for (size_t i = 0; i < y.dim(); i++)
        {
            for (const auto& e : A[i])
            {
                y[i] += alpha * e.val * x[e.idx];
            }
        }
    }
    else
    {
        for (size_t i = 0; i < x.dim(); i++)
        {
            for (const auto& e : A[i])
            {
                y[e.idx] += alpha * e.val * x[i];
            }
        }
    }
}
}
