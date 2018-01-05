#include "splinalg.hpp"
#include <cassert>
#include "blas.hpp"
#include "linalg.hpp"
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
real_t abs_max_val(vector_const_view v)
{
    size_t idx = blas_abs_max_idx(v);
    return abs(v[idx]);
}

real_t eigen_power_method(spmatrix_const_view A, vector_mutable_view x,
                          real_t tol, uint_t max_iter, uint_t check_interval)
{
    set_norm1(x, 1.0);
    vector x_next_(x.dim(), 0.0);
    vector_mutable_view x_next(x_next_);
    for (uint_t ii = 0; ii < max_iter / check_interval; ii++)
    {
        for (uint_t jj = 0; jj < check_interval; jj++)
        {
            dot(x_next, A, x);
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

real_t spsolve_sor_method(spmatrix_const_view A, vector_mutable_view x,
                          vector_const_view b, real_t w, real_t tol,
                          uint_t max_iter, uint_t check_interval)
{
    assert(A.m() == A.n());
    assert(A.m() == x.dim());
    assert(A.m() == b.dim());
    assert(A.is_compressed_row());
    vector x_next_(x.dim(), 0.0);
    vector_mutable_view x_next = x_next_;
    vector res(x.dim());
    real_t prec = 0;
    for (uint_t ii = 0; ii < max_iter / check_interval; ii++)
    {
        for (uint_t jj = 0; jj < check_interval; jj++)
        {
            for (size_t i = 0; i < A.ldim(); i++)
            {
                real_t diag = 0;
                real_t remain = b[i];
                auto view = A[i];
                for (size_t j = 0; j < view.nnz; j++)
                {
                    size_t idx = view.idx[j];
                    real_t val = view.val[j];
                    if (idx < i)
                    {
                        remain -= val * x_next[idx];
                    }
                    else if (idx > i)
                    {
                        remain -= val * x[idx];
                    }
                    else
                    {
                        diag = val;
                    }
                }
                assert(diag != 0);
                x_next[i] = (1 - w) * x[i] + w / diag * remain;
            }
            std::swap(x_next, x);
        }
        copy(res, b);
        spblas_matrix_vector(-1.0, A, x, 1.0, res);
        prec = abs_max_val(res);
        if (prec < tol)
        {
            return prec;
        }
    }
    return prec;
}

real_t spsolve_sor_method(spmatrix_const_view A, vector_mutable_view x,
                          vector_const_view b, real_t x_sum, real_t w,
                          real_t tol, uint_t max_iter, uint_t check_interval)
{
    assert(A.m() == A.n());
    assert(A.m() == x.dim());
    assert(A.m() == b.dim());
    assert(A.is_compressed_row());
    vector x_next_(x.dim(), 0.0);
    vector_mutable_view x_next = x_next_;
    size_t last_vec = A.ldim() - 1;
    vector res(x.dim());
    real_t prec = 0;
    for (uint_t ii = 0; ii < max_iter / check_interval; ii++)
    {
        for (uint_t jj = 0; jj < check_interval; jj++)
        {
            real_t last_remain = x_sum;
            for (size_t i = 0; i < last_vec; i++)
            {
                real_t diag = 0;
                real_t remain = b[i];
                auto view = A[i];
                for (size_t j = 0; j < view.nnz; j++)
                {
                    size_t idx = view.idx[j];
                    real_t val = view.val[j];
                    if (idx < i)
                    {
                        remain -= val * x_next[idx];
                    }
                    else if (idx > i)
                    {
                        remain -= val * x[idx];
                    }
                    else
                    {
                        diag = val;
                    }
                }
                x_next[i] = (1 - w) * x[i] + w / diag * remain;
                last_remain -= x_next[i];
            }
            x_next[last_vec] = (1 - w) * x[last_vec] + w * last_remain;
            std::swap(x_next, x);
        }
        copy(res, b);
        spblas_matrix_vector(-1.0, A, x, 1.0, res);
        prec = abs_max_val(res);
        if (prec < tol)
        {
            return prec;
        }
    }
    return prec;
}

void compute_xm(vector_mutable_view xm, vector_const_view x0,
                const std::vector<vector>& Vm, vector_const_view ym)
{
    assert(xm.dim() == x0.dim());
    assert(xm.dim() == Vm[0].dim());
    assert(Vm.size() == ym.dim());
    copy(xm, x0);
    for (size_t i = 0; i < ym.dim(); i++)
    {
        blas_axpy(ym[i], Vm[i], xm);
    }
}

// solve Ax = b, A must be full rank.
real_t spsolve_gmres_gms(spmatrix_const_view A, vector_mutable_view x,
                         vector_const_view b,
                         size_t kdim,  // dim of krylov space
                         real_t tol, uint_t check_interval)
{
    assert(A.m() == A.n());
    if (kdim > A.n())
    {
        kdim = A.n();
    }
    size_t M = A.m();
    std::vector<vector> V;
    V.reserve(kdim);
    matrix H(kdim + 1, kdim);
    // step 1a: V.col(0) = b - A * x
    V.push_back(vector(b));
    spblas_matrix_vector(-1.0, A, x, 1.0, V[0]);

    real_t beta = norm2(V[0]);
    scale(V[0], 1.0 / beta);
    vector wj(M);
    vector y(kdim + 1, 0.0);
    vector xm(x.dim());
    vector res(x.dim());
    real_t prec = 0;
    uint_t check_counter = 0;
    for (size_t j = 0; j < kdim; j++)
    {
        dot(wj, A, V[j]);
        for (size_t i = 0; i <= j; i++)
        {
            H(i, j) = dot(wj, V[i]);
            // wj = wj - H(i, j) * v_i
            blas_axpy(-H(i, j), V[i], wj);
        }
        H(j + 1, j) = norm2(wj);

        check_counter += 1;
        if ((check_counter % check_interval == 0) || (j == kdim - 1))
        {
            check_counter = 0;
            size_t m = j + 1;
            vector_mutable_view ym = y.sub(0, m + 1);
            fill(ym, 0);
            ym[0] = beta;
            auto Hm = matrix(H.sub(0, 0, m + 1, m));
            least_square_qr(Hm, ym);
            compute_xm(xm, x, V, ym.sub(0, m));
            copy(res, b);
            spblas_matrix_vector(-1.0, A, xm, 1.0, res);
            prec = abs_max_val(res);
            if (prec < tol || near_zero(H(j + 1, j), tol))
            {
                break;
            }
        }
        if (j < kdim - 1)
        {
            scale(wj, 1.0 / H(j + 1, j));
            V.push_back(wj);
        }
    }
    copy(x, xm);
    return prec;
}

real_t spsolve_restart_gmres_gms(spmatrix_const_view A, vector_mutable_view x,
                                 vector_const_view b,
                                 size_t kdim,  // dim of krylov space
                                 real_t tol, uint_t max_iter,
                                 uint_t check_interval)
{
    real_t prec = 0;
    for (uint_t ii = 0; ii < max_iter; ii++)
    {
        prec = spsolve_gmres_gms(A, x, b, kdim, tol, check_interval);
        if (prec < tol)
        {
            return prec;
        }
    }
    return prec;
}

real_t spdecomp_ilu(spmatrix_mutable_view A)
{
    assert(A.m() == A.n());
    assert(A.is_compressed_row());
    return 0;
}
}
