#include <cassert>
#include <sstream>
#include "blas.hpp"
#include "linalg.hpp"
#include "matvec_oper.hpp"
#include "spblas.hpp"
#include "splinalg.hpp"
#include "spmatvec_oper.hpp"
#include "vector.hpp"

namespace scigg
{
real_t abs_max_val(vector_const_view v)
{
    size_t idx = blas_abs_max_idx(v);
    return abs(v[idx]);
}

void spsolve_upper_tri(spmatrix_const_view U, vector_mutable_view b)
{
    assert(U.m() == U.n());
    assert(U.m() == b.dim());
    size_t N = U.n();
    std::vector<const real_t*> diag_ptr(N, nullptr);
    for (size_t i = 0; i < N; i++)
    {
        diag_ptr[i] = search_spvec_entry(U[i], i);
        if (diag_ptr[i] == nullptr || *diag_ptr[i] == 0.0)
        {
            throw matrix_value_error();
        }
    }
    if (U.is_compressed_row())
    {
        b[N - 1] /= *diag_ptr[N - 1];
        for (size_t ii = 1; ii < N; ii++)
        {
            size_t i = N - ii - 1;
            auto row = U[i];
            size_t skipped = static_cast<size_t>(diag_ptr[i] - row.val + 1);
            spvec_const_view subrow(row.nnz - skipped, row.idx + skipped,
                                    row.val + skipped);
            b[i] -= dot(subrow, b);
            b[i] /= *diag_ptr[i];
        }
    }
    else
    {
        b[N - 1] /= *diag_ptr[N - 1];
        for (size_t ii = 1; ii < N; ii++)
        {
            size_t i = N - ii - 1;
            auto col = U[i + 1];
            size_t count = static_cast<size_t>(diag_ptr[i + 1] - col.val);
            spvec_const_view subcol(count, col.idx, col.val);
            spblas_axpy(-b[i + 1], subcol, b);
            b[i] /= *diag_ptr[i];
        }
    }
}

void spsolve_lower_tri(spmatrix_const_view L, vector_mutable_view b)
{
    assert(L.m() == L.n());
    assert(L.m() == b.dim());
    size_t N = L.n();
    std::vector<const real_t*> diag_ptr(N, nullptr);
    for (size_t i = 0; i < N; i++)
    {
        diag_ptr[i] = search_spvec_entry(L[i], i);
        if (diag_ptr[i] == nullptr || *diag_ptr[i] == 0.0)
        {
            throw matrix_value_error();
        }
    }
    if (L.is_compressed_row())
    {
        b[0] /= *diag_ptr[0];
        for (size_t i = 1; i < N; i++)
        {
            auto row = L[i];
            size_t count = static_cast<size_t>(diag_ptr[i] - row.val);
            // entries before diagonal
            spvec_const_view subrow(count, row.idx, row.val);
            b[i] -= dot(subrow, b);
            b[i] /= *diag_ptr[i];
        }
    }
    else
    {
        b[0] /= *diag_ptr[0];
        for (size_t i = 1; i < N; i++)
        {
            auto col = L[i - 1];
            size_t skipped = static_cast<size_t>(diag_ptr[i - 1] - col.val + 1);
            spvec_const_view subcol(col.nnz - skipped, col.idx + skipped,
                                    col.val + skipped);
            spblas_axpy(-b[i - 1], subcol, b);
            b[i] /= *diag_ptr[i];
        }
    }
}

void spsolve_lower_tri_diag1(spmatrix_const_view L, vector_mutable_view b)
{
    assert(L.m() == L.n());
    assert(L.m() == b.dim());
    size_t N = L.n();
    std::vector<const size_t*> diag_idx_ptr(N, nullptr);
    for (size_t i = 0; i < N; i++)
    {
        auto vec = L[i];
        diag_idx_ptr[i] = std::lower_bound(vec.idx, vec.idx + vec.nnz, i);
    }
    if (L.is_compressed_row())
    {
        for (size_t i = 1; i < N; i++)
        {
            auto row = L[i];
            size_t count = static_cast<size_t>(diag_idx_ptr[i] - row.idx);
            // entries before diagonal
            spvec_const_view subrow(count, row.idx, row.val);
            b[i] -= dot(subrow, b);
        }
    }
    else
    {
        for (size_t i = 1; i < N; i++)
        {
            auto col = L[i - 1];
            size_t skipped =
                static_cast<size_t>(diag_idx_ptr[i - 1] - col.idx + 1);
            spvec_const_view subcol(col.nnz - skipped, col.idx + skipped,
                                    col.val + skipped);
            spblas_axpy(-b[i - 1], subcol, b);
        }
    }
}

void spsolve_ilu(spmatrix_const_view iLU, vector_mutable_view b)
{
    spsolve_lower_tri_diag1(iLU, b);
    spsolve_upper_tri(iLU, b);
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
    real_t prec = 0.0;
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
        copy(b, res);
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
    copy(x0, xm);
    for (size_t i = 0; i < ym.dim(); i++)
    {
        blas_axpy(ym[i], Vm[i], xm);
    }
}

// solve Ax = b, A must be full rank.
real_t spsolve_gmres_gms(spmatrix_const_view A, vector_mutable_view x,
                         vector_const_view b,
                         size_t kdim,  // dim of krylov space
                         real_t tol, uint_t check_interval,
                         pre_condition Msolve)
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
    // step 1a: V.col(0) = M^-1 (b - A * x)
    V.push_back(vector(b));
    spblas_matrix_vector(-1.0, A, x, 1.0, V[0]);
    if (Msolve)
    {
        Msolve(V[0]);
    }

    real_t beta = norm2(V[0]);
    scale(1.0 / beta, V[0]);
    vector wj(M);
    vector y(kdim + 1, 0.0);
    vector xm(x.dim());
    vector res(x.dim());
    real_t prec = 0;
    uint_t check_counter = 0;
    for (size_t j = 0; j < kdim; j++)
    {
        // wj = M^-1 A * V[j]
        dot(A, V[j], wj);
        if (Msolve)
        {
            Msolve(wj);
        }
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
            fill(0, ym);
            ym[0] = beta;
            auto Hm = matrix(H.sub(0, 0, m + 1, m));
            least_square_qr(Hm, ym);
            compute_xm(xm, x, V, ym.sub(0, m));
            copy(b, res);
            // res = (b - A xm), no precondition for res
            spblas_matrix_vector(-1.0, A, xm, 1.0, res);
            prec = abs_max_val(res);
            if (prec < tol || near_zero(H(j + 1, j), tol))
            {
                break;
            }
        }
        if (j < kdim - 1)
        {
            scale(1.0 / H(j + 1, j), wj);
            V.push_back(wj);
        }
    }
    copy(xm, x);
    return prec;
}

real_t spsolve_restart_gmres_gms(spmatrix_const_view A, vector_mutable_view x,
                                 vector_const_view b,
                                 size_t kdim,  // dim of krylov space
                                 real_t tol, uint_t max_iter,
                                 uint_t check_interval, pre_condition Msolve)
{
    real_t prec = 0;
    for (uint_t ii = 0; ii < max_iter; ii++)
    {
        prec = spsolve_gmres_gms(A, x, b, kdim, tol, check_interval, Msolve);
        if (prec < tol)
        {
            return prec;
        }
    }
    return prec;
}
}
