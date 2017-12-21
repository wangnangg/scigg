#include "matvec_oper.hpp"
#include "blas.hpp"
namespace markovgg
{
matrix create_matrix(size_t m, size_t n, const std::vector<double>& v)
{
    assert(m * n == v.size());
    matrix M(m, n);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            M(i, j) = v[i * n + j];
        }
    }
    return M;
}
vector create_vector(size_t n, const std::vector<double>& v)
{
    vector vo(n);
    for (size_t i = 0; i < n; i++)
    {
        vo[i] = v[i];
    }
    return vo;
}
vector_const_view row_const_view(matrix_const_view A, size_t row_idx)
{
    return vector_const_view(&A(row_idx, 0), A.n(), 1);
}

vector_mutable_view row_mutable_view(matrix_mutable_view A, size_t row_idx)
{
    return vector_mutable_view(&A(row_idx, 0), A.n(), 1);
}

vector_const_view col_const_view(matrix_const_view A, size_t col_idx)
{
    return vector_const_view(&A(0, col_idx), A.m(), A.ldim());
}

vector_mutable_view col_mutable_view(matrix_mutable_view A, size_t col_idx)
{
    return vector_mutable_view(&A(0, col_idx), A.m(), A.ldim());
}

bool near_zero(vector_const_view v1, real_t tol)
{
    size_t idx = blas_abs_max_idx(v1);
    if (near_zero(v1[idx], tol))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool near_eq(vector_const_view v1, vector_const_view v2, real_t tol)
{
    assert(v1.dim() == v2.dim());
    real_t max_v = abs(v1[0] - v2[0]);
    for (size_t i = 1; i < v1.dim(); i++)
    {
        real_t v = abs(v1[i] - v2[i]);
        if (max_v < v)
        {
            max_v = v;
        }
    }
    return near_zero(max_v, tol);
}

bool near_eq(matrix_const_view M1, matrix_const_view M2, real_t tol)
{
    assert(M1.m() == M2.m());
    assert(M1.n() == M2.n());
    for (size_t i = 0; i < M1.m(); i++)
    {
        for (size_t j = 0; j < M1.n(); j++)
        {
            if (!near_zero(M1(i, j) - M2(i, j), tol))
            {
                return false;
            }
        }
    }
    return true;
}

// C = A * B
void dot(matrix_mutable_view C, matrix_const_view A, bool transposeA,
         matrix_const_view B, bool transposeB)
{
    blas_matrix_matrix(1.0, A, transposeA, B, transposeB, 0.0, C);
}
matrix dot(matrix_const_view A, bool transposeA, matrix_const_view B,
           bool transposeB)
{
    auto C = matrix(A.m(), B.n());
    blas_matrix_matrix(1.0, A, transposeA, B, transposeB, 0.0, C);
    return C;
}
}
