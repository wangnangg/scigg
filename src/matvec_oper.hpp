#pragma once
#include "matrix.hpp"
#include "vector.hpp"

namespace markovgg
{
matrix create_matrix(size_t m, size_t n, const std::vector<double>& v);
vector create_vector(size_t n, const std::vector<double>& v);

vector_const_view matrix_row_const_view(matrix_const_view A, size_t row_idx,
                                        size_t start_col = 0,
                                        size_t end_col = 0);
vector_mutable_view matrix_row_mutable_view(matrix_mutable_view A,
                                            size_t row_idx,
                                            size_t start_col = 0,
                                            size_t end_col = 0);
vector_const_view matrix_col_const_view(matrix_const_view A, size_t col_idx,
                                        size_t start_row = 0,
                                        size_t end_row = 0);
vector_mutable_view matrix_col_mutable_view(matrix_mutable_view A,
                                            size_t col_idx,
                                            size_t start_row = 0,
                                            size_t end_row = 0);

matrix_mutable_view submatrix_mutable_view(matrix_mutable_view A,
                                           size_t start_row, size_t start_col,
                                           size_t end_row = 0,
                                           size_t end_col = 0);
matrix_const_view submatrix_const_view(matrix_const_view A, size_t start_row,
                                       size_t start_col, size_t end_row = 0,
                                       size_t end_col = 0);
vector_const_view subvector_const_view(vector_const_view v, size_t start_idx,
                                       size_t end_idx = 0);
vector_mutable_view subvector_mutable_view(vector_mutable_view v,
                                           size_t start_idx,
                                           size_t end_idx = 0);

inline real_t abs(real_t v) { return v >= 0 ? v : -v; }
inline real_t max(size_t x, size_t y) { return x > y ? x : y; }
inline real_t min(size_t x, size_t y) { return x > y ? y : x; }

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

void fill(vector_mutable_view x, real_t val);

void set_norm1(vector_mutable_view vec, real_t val);
void set_norm2(vector_mutable_view vec, real_t val);

// y = x + y
void inc(vector_mutable_view y, vector_const_view x);
// y = x - y
void dec(vector_mutable_view y, vector_const_view x);
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
