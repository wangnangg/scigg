#include "linalg.hpp"
#include <cassert>
#include "blas.hpp"
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
        auto qi = matrix_col_mutable_view(A, i);
        real_t rii = blas_norm2(qi);
        blas_scale(1.0 / rii, qi);
        R(i, i) = rii;
        for (size_t j = i + 1; j < A.n(); j++)
        {
            auto vj = matrix_col_mutable_view(A, j);
            real_t rij = blas_dot(qi, vj);
            // vj = vj - rjj * qi
            blas_axpy(-rij, qi, vj);
            R(i, j) = rij;
        }
    }
}

void qr_decomp_hr(matrix_mutable_view A, matrix_mutable_view V)
{
    assert(A.m() >= A.n());
    assert(A.m() == V.m());
    assert(A.n() == V.n());
    size_t m = A.m();
    size_t n = A.n();
    vector vks_A_(n);  // workspace
    for (size_t k = 0; k < n; k++)
    {
        vector_const_view xk = matrix_col_const_view(A, k, k, m);
        vector_mutable_view vk = matrix_col_mutable_view(V, k, k, m);
        blas_copy(xk, vk);
        real_t xnorm = blas_norm2(xk);
        if (xk[0] < 0)
        {
            xnorm = -xnorm;
        }
        vk[0] += xnorm;
        set_norm2(vk, 1.0);
        vector_mutable_view vks_Ak = subvector_mutable_view(vks_A_, k);
        matrix_mutable_view Ak = submatrix_mutable_view(A, k, k);
        // vks_Ak = vk' . Ak
        blas_matrix_vector(1.0, Ak, true, vk, 0.0, vks_Ak);
        // rank 1 update: Ak = -2 * vk * vks_Ak  + Ak
        blas_rank1(-2.0, vk, vks_Ak, Ak);
    }
}

void recover_q_from_v(matrix_const_view V, matrix_mutable_view Q)
{
    assert(V.m() == Q.m());
    assert(Q.m() == Q.n());
    size_t m = V.m();
    size_t n = V.n();
    for (size_t i = 0; i < m; i++)
    {
        vector_mutable_view qi = matrix_col_mutable_view(Q, i);
        fill(qi, 0.0);
        qi[i] = 1.0;
        size_t num_ref = min(n, i + 1);
        for (size_t j = 0; j < num_ref; j++)
        {
            size_t k = num_ref - j - 1;
            auto vk = matrix_col_const_view(V, k, k, m);
            auto qi_k = subvector_mutable_view(qi, k, m);
            real_t alpha = -2 * blas_dot(qi_k, vk);
            blas_axpy(alpha, vk, qi_k);
        }
    }
}
}
