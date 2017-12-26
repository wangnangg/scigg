#include "splinalg.hpp"
#include <cassert>
#include "blas.hpp"
#include "matvec_oper.hpp"
#include "spblas.hpp"
#include "spmatrix_oper.hpp"
#include "vector.hpp"

namespace markovgg
{
real_t max_diff(vector_const_view v1, vector_const_view v2)
{
    assert(v1.dim() == v2.dim());
    real_t max = 0;
    for (size_t i = 0; i < v1.dim(); i++)
    {
        real_t diff = abs(v1[i] - v2[i]);
        if (diff > max)
        {
            max = diff;
        }
    }
    return max;
}

real_t eigen_power_method(const spmatrix& A, bool transposeA,
                          vector_mutable_view x, real_t tol, int_t max_iter,
                          int_t check_interval)
{
    set_norm1(x, 1.0);
    vector x_next_(x.dim(), 0.0);
    vector_mutable_view x_next(x_next_);
    for (int_t ii = 0; ii < max_iter / check_interval; ii++)
    {
        for (int_t jj = 0; jj < check_interval; jj++)
        {
            dot(x_next, A, transposeA, x);
            set_norm1(x_next, 1.0);
            std::swap(x_next, x);
        }
        real_t prec = max_diff(x_next, x);
        if (prec < tol)
        {
            return prec;
        }
    }
    return max_diff(x_next, x);
}

real_t linsv_sor_method(const spmatrix& A, bool transposeA,
                        vector_mutable_view x, vector_const_view b, real_t w,
                        real_t tol, int_t max_iter, int_t check_interval)
{
    assert(A.m() == A.n());
    assert(A.m() == x.dim());
    assert(A.m() == b.dim());
    assert((!transposeA && A.format() == CPR_ROW) ||
           (transposeA && A.format() == CPR_COL));
    vector x_next_(x.dim(), 0.0);
    vector_mutable_view x_next = x_next_;
    for (int_t ii = 0; ii < max_iter / check_interval; ii++)
    {
        for (int_t jj = 0; jj < check_interval; jj++)
        {
            for (size_t i = 0; i < A.ldim(); i++)
            {
                real_t diag = 0;
                real_t remain = b[i];
                for (const auto& e : A[i])
                {
                    if (e.idx < i)
                    {
                        remain -= e.val * x_next[e.idx];
                    }
                    else if (e.idx > i)
                    {
                        remain -= e.val * x[e.idx];
                    }
                    else
                    {
                        diag = e.val;
                    }
                }
                x_next[i] = (1 - w) * x[i] + w / diag * remain;
            }
            std::swap(x_next, x);
        }
        real_t prec = max_diff(x_next, x);
        if (prec < tol)
        {
            return prec;
        }
    }
    return max_diff(x_next, x);
}

real_t linsv_sor_method(const spmatrix& A, bool transposeA,
                        vector_mutable_view x, vector_const_view b,
                        real_t x_sum, real_t w, real_t tol, int_t max_iter,
                        int_t check_interval)
{
    assert(A.m() == A.n());
    assert(A.m() == x.dim());
    assert(A.m() == b.dim());
    assert((!transposeA && A.format() == CPR_ROW) ||
           (transposeA && A.format() == CPR_COL));
    vector x_next_(x.dim(), 0.0);
    vector_mutable_view x_next = x_next_;
    size_t last_vec = A.ldim() - 1;
    for (int_t ii = 0; ii < max_iter / check_interval; ii++)
    {
        for (int_t jj = 0; jj < check_interval; jj++)
        {
            real_t last_remain = x_sum;
            for (size_t i = 0; i < last_vec; i++)
            {
                real_t diag = 0;
                real_t remain = b[i];
                for (const auto& e : A[i])
                {
                    if (e.idx < i)
                    {
                        remain -= e.val * x_next[e.idx];
                    }
                    else if (e.idx > i)
                    {
                        remain -= e.val * x[e.idx];
                    }
                    else
                    {
                        diag = e.val;
                    }
                }
                x_next[i] = (1 - w) * x[i] + w / diag * remain;
                last_remain -= x_next[i];
            }
            x_next[last_vec] = (1 - w) * x[last_vec] + w * last_remain;
            std::swap(x_next, x);
        }
        real_t prec = max_diff(x_next, x);
        if (prec < tol)
        {
            return prec;
        }
    }
    return max_diff(x_next, x);
}

// solve Ax = b, A must be full rank.
real_t linsv_gmres(const spmatrix& A, bool transposeA, vector_mutable_view x,
                   vector_const_view b,
                   size_t m,  // dim of krylov space
                   real_t tol)
{
    assert(A.m() == A.n());
    size_t n = A.n();
    matrix H(m, n);
    return 0;
}
}
