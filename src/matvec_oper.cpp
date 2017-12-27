#include "matvec_oper.hpp"
#include "blas.hpp"
namespace markovgg
{
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

void fill(vector_mutable_view x, real_t val)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] = val;
    }
}

void fill(matrix_mutable_view A, real_t val)
{
    for (size_t i = 0; i < A.m(); i++)
    {
        for (size_t j = 0; j < A.n(); j++)
        {
            A(i, j) = val;
        }
    }
}

void set_diag(matrix_mutable_view A, real_t val)
{
    fill(A, 0.0);
    size_t iter = min(A.m(), A.n());
    for (size_t i = 0; i < iter; i++)
    {
        A(i, i) = 1.0;
    }
}

vector operator*(real_t alpha, vector_const_view x)
{
    vector y(x.dim(), 0.0);
    blas_axpy(alpha, x, y);
    return y;
}

void scale(matrix_mutable_view A, real_t alpha)
{
    for (size_t i = 0; i < A.m(); i++)
    {
        for (size_t j = 0; j < A.n(); j++)
        {
            A(i, j) *= alpha;
        }
    }
}

matrix operator*(real_t alpha, matrix_const_view A)
{
    matrix B(A.m(), A.n());
    for (size_t i = 0; i < B.m(); i++)
    {
        for (size_t j = 0; j < B.n(); j++)
        {
            B(i, j) = A(i, j) * alpha;
        }
    }
    return B;
}
void add(vector_mutable_view z, vector_const_view x, vector_const_view y)
{
    assert(z.dim() == x.dim());
    assert(z.dim() == y.dim());
    for (size_t i = 0; i < z.dim(); i++)
    {
        z[i] = x[i] + y[i];
    }
}

void inc(vector_mutable_view z, vector_const_view x)
{
    // z = x + z
    blas_axpy(1.0, x, z);
}

void sub(vector_mutable_view z, vector_const_view x, vector_const_view y)
{
    assert(z.dim() == x.dim());
    assert(z.dim() == y.dim());
    for (size_t i = 0; i < z.dim(); i++)
    {
        z[i] = x[i] - y[i];
    }
}

void dec(vector_mutable_view z, vector_const_view x)
{
    // z = - x + z
    blas_axpy(-1.0, x, z);
}
real_t dot(vector_const_view x, vector_const_view y) { return blas_dot(x, y); }

void dot(vector_mutable_view z, matrix_const_view A, vector_const_view x)
{
    blas_matrix_vector(1.0, A, x, 0.0, z);
}
// C = A * B
void dot(matrix_mutable_view C, matrix_const_view A, matrix_const_view B)
{
    blas_matrix_matrix(1.0, A, B, 0.0, C);
}
}
