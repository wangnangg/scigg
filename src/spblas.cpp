#include "spblas.hpp"
#include "matvec_oper.hpp"
namespace markovgg
{
// y = alpha * A * x + beta * y
void spblas_matrix_vector(real_t alpha, spmatrix_const_view A,
                          vector_const_view x, real_t beta,
                          vector_mutable_view y)
{
    assert(A.n() == x.dim());
    assert(A.m() == y.dim());
    scale(y, beta);
    if (A.is_compressed_row())
    {
        for (size_t i = 0; i < y.dim(); i++)
        {
            auto view = A[i];
            for (size_t j = 0; j < view.nnz; j++)
            {
                y[i] += alpha * view.val[j] * x[view.idx[j]];
            }
        }
    }
    else
    {
        for (size_t i = 0; i < x.dim(); i++)
        {
            auto view = A[i];
            for (size_t j = 0; j < view.nnz; j++)
            {
                y[view.idx[j]] += alpha * view.val[j] * x[i];
            }
        }
    }
}
}
