#include <cassert>
#include "blas.hpp"
#include "debug_utils.hpp"
#include "linalg.hpp"
#include "matvec_oper.hpp"

namespace markovgg
{
void decomp_qr_mgs(matrix_mutable_view A, matrix_mutable_view R)
{
    assert(A.m() >= A.n());
    assert(R.m() == R.n());
    assert(A.n() == R.m());
    for (size_t i = 0; i < A.n(); i++)
    {
        auto qi = A.col(i);
        real_t rii = norm2(qi);
        scale(1.0 / rii, qi);
        R(i, i) = rii;
        for (size_t j = i + 1; j < A.n(); j++)
        {
            auto vj = A.col(j);
            real_t rij = dot(qi, vj);
            // vj = vj - rjj * qi
            blas_axpy(-rij, qi, vj);
            R(i, j) = rij;
        }
    }
}

// find vector v, with v[0] = 1, so that (I - tau v v') w zeros w[1:n]. v[1:n]
// is stored in w[1:n] and w[0] = (Px)[0]
real_t find_householder_vector(vector_mutable_view w)
{
    real_t wr = norm2(w);
    if (w[0] < 0)
    {
        wr = -wr;
    }
    w[0] += wr;
    scale(1.0 / w[0], w);
    real_t vnorm = norm2(w);
    real_t tau = 2.0 / (vnorm * vnorm);
    w[0] = -wr;
    return tau;
}

// w = w - \tau vv'w, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(vector_mutable_view w, real_t tau,
                                 vector_const_view v)
{
    auto vm = vector_mutable_view(const_cast<real_t*>(&v[0]), v.dim(), v.inc());
    real_t v0 = vm[0];
    vm[0] = 1;
    real_t val = dot(vm, w);
    blas_axpy(-tau * val, vm, w);
    vm[0] = v0;
}

// A = A - \tau vv'A, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(matrix_mutable_view A, real_t tau,
                                 vector_const_view v)
{
    auto vm = vector_mutable_view(const_cast<real_t*>(&v[0]), v.dim(), v.inc());
    real_t v0 = vm[0];
    vm[0] = 1;
    auto vt_A = vm * A;
    blas_rank1(-tau, vm, vt_A, A);
    vm[0] = v0;
}

void decomp_qr_hr(matrix_mutable_view A, vector_mutable_view tau_vec)
{
    assert(A.m() >= A.n());
    assert(tau_vec.dim() == A.n());
    size_t n = A.n();
    for (size_t k = 0; k < n; k++)
    {
        auto wk = A.col(k).sub(k);
        real_t tau = find_householder_vector(wk);
        if (k < n - 1)
        {
            auto Ak = A.sub(k, k + 1);
            apply_householder_reflector(Ak, tau, wk);
        }
        tau_vec[k] = tau;
    }
}

void unpack_qr(matrix_const_view QR, vector_const_view tau_vec,
               matrix_mutable_view Q, matrix_mutable_view R)
{
    assert(QR.m() == Q.m());
    assert(QR.n() == Q.n());
    assert(R.m() == R.n());
    assert(R.m() == QR.n());
    assert(QR.n() == tau_vec.dim());
    size_t n = QR.n();
    fill(0.0, Q);
    for (size_t i = 0; i < Q.n(); i++)
    {
        Q(i, i) = 1.0;
    }
    for (size_t i = 0; i < n; i++)
    {
        size_t k = n - i - 1;
        auto vk = QR.col(k).sub(k);
        auto Qk = Q.sub(k, k);
        apply_householder_reflector(Qk, tau_vec[k], vk);
    }
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            R(i, j) = 0.0;
        }
        for (size_t j = i; j < n; j++)
        {
            R(i, j) = QR(i, j);
        }
    }
}

// compute v = Q^T * v
void qt_dot_vector(matrix_const_view QR, vector_const_view tau,
                   vector_mutable_view v)
{
    assert(QR.m() >= QR.n());
    assert(QR.n() == tau.dim());
    assert(QR.m() == v.dim());
    for (size_t i = 0; i < QR.n(); i++)
    {
        apply_householder_reflector(v.sub(i), tau[i], QR.col(i).sub(i));
    }
}

// compute v = Q * v
void q_dot_vector(matrix_const_view QR, vector_const_view tau,
                  vector_mutable_view v)
{
    assert(QR.m() >= QR.n());
    assert(QR.n() == tau.dim());
    assert(QR.m() == v.dim());
    for (size_t i = 0; i < QR.n(); i++)
    {
        size_t k = QR.n() - i - 1;
        apply_householder_reflector(v.sub(k), tau[k], QR.col(k).sub(k));
    }
}

// PA = LU, A will be replaced by L and U, the rows of P will be permuted
void decomp_lup(matrix_mutable_view A, matrix_mutable_view P)
{
    assert(A.m() == A.n());
    assert(P.m() == A.m());
    size_t N = A.n();
    for (size_t j = 0; j < N - 1; j++)
    {
        size_t pivot_idx = blas_abs_max_idx(A.col(j).sub(j)) + j;
        if (j != pivot_idx)
        {
            blas_swap(A.row(j), A.row(pivot_idx));
            blas_swap(P.row(j), P.row(pivot_idx));
        }
        auto rowj = A.row(j);
        real_t pivot_val = A(j, j);
        if (pivot_val == 0.0)
        {
            throw matrix_value_error();
        }
        for (size_t i = j + 1; i < N; i++)
        {
            real_t ratio = A(i, j) / pivot_val;
            if (ratio != 0)
            {
                A(i, j) = ratio;
                blas_axpy(-ratio, rowj.sub(j + 1), A.row(i).sub(j + 1));
            }
        }
    }
}

void decomp_lu(matrix_mutable_view A)
{
    assert(A.m() == A.n());
    size_t N = A.n();
    for (size_t j = 0; j < N - 1; j++)
    {
        auto rowj = A.row(j);
        real_t pivot_val = A(j, j);
        if (pivot_val == 0.0)
        {
            throw matrix_value_error();
        }
        for (size_t i = j + 1; i < N; i++)
        {
            real_t ratio = A(i, j) / pivot_val;
            if (ratio != 0)
            {
                A(i, j) = ratio;
                blas_axpy(-ratio, rowj.sub(j + 1), A.row(i).sub(j + 1));
            }
        }
    }
}

// A will be replaced by U
void unpack_lu(matrix_mutable_view A, matrix_mutable_view L)
{
    assert(A.m() == A.n());
    assert(L.m() == L.n());
    assert(A.m() == L.m());
    size_t M = A.m();
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            L(i, j) = A(i, j);
            A(i, j) = 0.0;
        }
        L(i, i) = 1.0;
    }
}
}
