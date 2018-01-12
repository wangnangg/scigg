#include "spblas.hpp"
#include "matvec_oper.hpp"
namespace markovgg
{
void spblas_zero_fillin_axpy(real_t alpha, spvec_const_view x,
                             spvec_mutable_view y)
{
    size_t x_i = 0;
    size_t y_i = 0;
    while (x_i < x.nnz && y_i < y.nnz)
    {
        size_t xcol = x.idx[x_i];
        size_t ycol = y.idx[y_i];
        if (xcol == ycol)
        {
            y.val[y_i] += alpha * x.val[x_i];
            x_i += 1;
            y_i += 1;
        }
        else if (xcol < ycol)
        {
            x_i += 1;
        }
        else if (xcol > ycol)
        {
            y_i += 1;
        }
    }
}

void spblas_axpy(real_t alpha, spvec_const_view x, vector_mutable_view y)
{
    for (size_t i = 0; i < x.nnz; i++)
    {
        size_t idx = x.idx[i];
        real_t val = x.val[i];
        y[idx] += alpha * val;
    }
}

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
