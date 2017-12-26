#include "linalg.hpp"
#include <cassert>
#include "blas.hpp"
#include "debug_utils.hpp"
#include "matvec_oper.hpp"

namespace markovgg
{
void qr_decomp_mgs(matrix_mutable_view A, matrix_mutable_view R)
{
    assert(A.m() >= A.n());
    assert(R.m() == R.n());
    assert(A.n() == R.m());
    for (size_t i = 0; i < A.n(); i++)
    {
        auto qi = A.col(i);
        real_t rii = blas_norm2(qi);
        blas_scale(1.0 / rii, qi);
        R(i, i) = rii;
        for (size_t j = i + 1; j < A.n(); j++)
        {
            auto vj = A.col(j);
            real_t rij = blas_dot(qi, vj);
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
    real_t wr = blas_norm2(w);
    if (w[0] < 0)
    {
        wr = -wr;
    }
    w[0] += wr;
    blas_scale(1.0 / w[0], w);
    real_t vnorm = blas_norm2(w);
    real_t tau = 2.0 / (vnorm * vnorm);
    w[0] = -wr;
    return tau;
}

// w = w - \tau vv'w, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(vector_mutable_view w, real_t tau,
                                 vector_mutable_view v)
{
    real_t v0 = v[0];
    v[0] = 1;
    real_t val = blas_dot(v, w);
    blas_axpy(-tau * val, v, w);
    v[0] = v0;
}

// A = A - \tau vv'A, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(matrix_mutable_view A, real_t tau,
                                 vector_mutable_view v)
{
    real_t v0 = v[0];
    v[0] = 1;
    vector vt_A(A.n());
    blas_matrix_vector(1.0, A.transpose(), v, 0.0, vt_A);
    blas_rank1(-tau, v, vt_A, A);
    v[0] = v0;
}

void qr_decomp_hr(matrix_mutable_view A, vector_mutable_view tau_vec)
{
    assert(A.m() >= A.n());
    assert(tau_vec.dim() == A.n());
    size_t m = A.m();
    size_t n = A.n();
    vector vks_A_(n);  // workspace
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
void unpack_qr(matrix_mutable_view QR, vector_const_view tau_vec,
               matrix_mutable_view Q, matrix_mutable_view R)
{
    assert(QR.m() == Q.m());
    assert(QR.n() == Q.n());
    assert(R.m() == R.n());
    assert(R.m() == QR.n());
    assert(QR.n() == tau_vec.dim());
    size_t m = QR.m();
    size_t n = QR.n();
    fill(Q, 0.0);
    set_diag(Q, 1.0);
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
}
