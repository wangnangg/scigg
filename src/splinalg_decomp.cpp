#include <cassert>
#include <sstream>
#include "blas.hpp"
#include "linalg.hpp"
#include "matvec_oper.hpp"
#include "spblas.hpp"
#include "splinalg.hpp"
#include "spmatrix_oper.hpp"
#include "vector.hpp"

namespace markovgg
{
void zero_fillin_axpy(real_t alpha, spvec_const_view x, spvec_mutable_view y)
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

void spdecomp_ilu(spmatrix_mutable_view A)
{
    assert(A.m() == A.n());
    assert(A.is_compressed_row());
    size_t N = A.n();
    std::vector<real_t*> diag_ptr(N, nullptr);
    for (size_t i = 0; i < N; i++)
    {
        diag_ptr[i] = search_spvec_entry(A[i], i);
        if (diag_ptr[i] == nullptr)
        {
            throw matrix_value_error();
        }
    }
    if (*diag_ptr[0] == 0.0)
    {
        throw matrix_value_error();
    }
    for (size_t i = 1; i < N; i++)
    {
        auto row = A[i];
        // entries before diagonal
        for (size_t j = 0; (j < row.nnz) && (row.val + j < diag_ptr[i]); j++)
        {
            size_t idx = row.idx[j];
            real_t ratio = row.val[j] / (*diag_ptr[idx]);
            row.val[j] = ratio;
            zero_fillin_axpy(-ratio, A[idx].sub(idx + 1), row.sub(idx + 1));
        }
        if (*diag_ptr[i] == 0.0)
        {
            throw matrix_value_error();
        }
    }
}
}
