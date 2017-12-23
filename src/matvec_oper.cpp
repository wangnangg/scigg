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
vector_const_view matrix_row_const_view(matrix_const_view A, size_t row_idx,
                                        size_t start_col, size_t end_col)
{
    assert(end_col <= A.n());
    if (end_col == 0)
    {
        end_col = A.n();
    }
    size_t dim = end_col - start_col;
    return vector_const_view(&A(row_idx, start_col), dim, 1);
}

vector_mutable_view matrix_row_mutable_view(matrix_mutable_view A,
                                            size_t row_idx, size_t start_col,
                                            size_t end_col)
{
    assert(end_col <= A.n());
    if (end_col == 0)
    {
        end_col = A.n();
    }
    size_t dim = end_col - start_col;
    return vector_mutable_view(&A(row_idx, start_col), dim, 1);
}

vector_const_view matrix_col_const_view(matrix_const_view A, size_t col_idx,
                                        size_t start_row, size_t end_row)
{
    assert(end_row <= A.m());
    if (end_row == 0)
    {
        end_row = A.m();
    }
    size_t dim = end_row - start_row;
    return vector_const_view(&A(start_row, col_idx), dim, A.ldim());
}

vector_mutable_view matrix_col_mutable_view(matrix_mutable_view A,
                                            size_t col_idx, size_t start_row,
                                            size_t end_row)
{
    assert(end_row <= A.m());
    if (end_row == 0)
    {
        end_row = A.m();
    }
    size_t dim = end_row - start_row;
    return vector_mutable_view(&A(start_row, col_idx), dim, A.ldim());
}

matrix_mutable_view submatrix_mutable_view(matrix_mutable_view A,
                                           size_t start_row, size_t start_col,
                                           size_t end_row, size_t end_col)
{
    assert(end_row <= A.m());
    assert(end_col <= A.n());
    if (end_row == 0)
    {
        end_row = A.m();
    }
    if (end_col == 0)
    {
        end_col = A.n();
    }
    size_t m = end_row - start_row;
    size_t n = end_col - start_col;
    return matrix_mutable_view(&A(start_row, start_col), m, n, A.ldim());
}

matrix_const_view submatrix_const_view(matrix_const_view A, size_t start_row,
                                       size_t start_col, size_t end_row,
                                       size_t end_col)
{
    assert(end_row <= A.m());
    assert(end_col <= A.n());
    if (end_row == 0)
    {
        end_row = A.m();
    }
    if (end_col == 0)
    {
        end_col = A.n();
    }
    size_t m = end_row - start_row;
    size_t n = end_col - start_col;
    return matrix_const_view(&A(start_row, start_col), m, n, A.ldim());
}

vector_const_view subvector_const_view(vector_const_view v, size_t start_idx,
                                       size_t end_idx)
{
    assert(end_idx <= v.dim());
    if (end_idx == 0)
    {
        end_idx = v.dim();
    }
    return vector_const_view(&v[start_idx], end_idx - start_idx, v.inc());
}
vector_mutable_view subvector_mutable_view(vector_mutable_view v,
                                           size_t start_idx, size_t end_idx)
{
    assert(end_idx <= v.dim());
    if (end_idx == 0)
    {
        end_idx = v.dim();
    }
    return vector_mutable_view(&v[start_idx], end_idx - start_idx, v.inc());
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

void inc(vector_mutable_view y, vector_const_view x) { blas_axpy(1.0, x, y); }
void dec(vector_mutable_view y, vector_const_view x) { blas_axpy(-1.0, x, y); }
void fill(vector_mutable_view x, real_t val)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] = val;
    }
}

void set_norm1(vector_mutable_view vec, real_t val)
{
    real_t n1 = blas_abs_sum(vec);
    blas_scale(val / n1, vec);
}

void set_norm2(vector_mutable_view vec, real_t val)
{
    real_t n2 = blas_norm2(vec);
    blas_scale(val / n2, vec);
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
