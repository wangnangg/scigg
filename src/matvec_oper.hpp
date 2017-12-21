#pragma once
#include "matrix.hpp"
#include "vector.hpp"

namespace markovgg
{
matrix create_matrix(size_t m, size_t n, const std::vector<double>& v);
vector create_vector(size_t n, const std::vector<double>& v);

vector_const_view row_const_view(matrix_const_view A, size_t row_idx);
vector_mutable_view row_mutable_view(matrix_mutable_view A, size_t row_idx);
vector_const_view col_const_view(matrix_const_view A, size_t col_idx);
vector_mutable_view col_mutable_view(matrix_mutable_view A, size_t col_idx);

template <typename mat_type>
vector row_vector(const mat_type& A, size_t row_idx)
{
    vector v(A.n());
    for (size_t i = 0; i < A.n(); i++)
    {
        v[i] = A(row_idx, i);
    }
    return v;
}

template <typename mat_type>
vector col_vector(const mat_type& A, size_t col_idx)
{
    vector v(A.m());
    for (size_t i = 0; i < A.m(); i++)
    {
        v[i] = A(i, col_idx);
    }
    return v;
}

inline real_t abs(real_t v)
{
    return v >= 0 ? v : -v;
    ;
}

inline bool near_zero(real_t v1, real_t tol)
{
    if (abs(v1) < tol)
    {
        return true;
    }
    else
    {
        return false;
    }
}
inline bool near_eq(real_t v1, real_t v2, real_t tol)
{
    return near_zero(v1 - v2, tol);
}

bool near_zero(const vector& v1, real_t tol);
bool near_eq(vector_const_view v1, vector_const_view v2, real_t tol);
bool near_eq(matrix_const_view M1, matrix_const_view M2, real_t tol);

// y = A . x
void dot(vector_mutable_view y, matrix_const_view A, bool transposeA,
         vector_mutable_view x);
vector dot(matrix_const_view A, bool transposeA, vector_mutable_view x);

// C = A . B
void dot(matrix_mutable_view C, matrix_const_view A, bool transposeA,
         matrix_const_view B, bool transposeB);
matrix dot(matrix_const_view A, bool transposeA, matrix_const_view B,
           bool transposeB);
}
