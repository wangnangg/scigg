#pragma once
#include <cmath>
#include "matrix.hpp"
namespace markovgg
{
inline bool near_zero(real_t v, real_t prec) { return absolute(v) <= prec; }
inline bool near_eq(real_t a, real_t b, real_t prec)
{
    return near_zero(a - b, prec);
}
bool near_eq(const vec& x, const vec& y, real_t prec);
bool near_zero(const vec& x, real_t prec);
real_t norm1(const vec& vec);
real_t norm2(const vec& vec);
void set_norm1(vec& vec, real_t val);
void set_norm2(vec& vec, real_t val);
real_t vec_sum(const vec& v);
void set_vec_sum(vec& v, real_t val);

void fill(vec& x, real_t val);

// x += v - w
void acc_sub(vec& x, const vec& v, const vec& w);
// x = v - w
void sub(vec& x, const vec& v, const vec& w);
inline vec sub(const vec& v, const vec& w)
{
    vec x(v.dim(), 0.0);
    sub(x, v, w);
    return x;
}
// x -= v
void self_add(vec& x, const vec& v);

// x += v + w
void acc_add(vec& x, const vec& v, const vec& w);

// x = v + w
void add(vec& x, const vec& v, const vec& w);
inline vec add(const vec& v, const vec& w)
{
    vec x(v.dim(), 0.0);
    add(x, v, w);
    return x;
}

// x += v . A
void acc_mul(vec& x, const vec& v, const matrix& A);
// x = v . A
void mul(vec& x, const vec& v, const matrix& A);
inline vec mul(const vec& v, const matrix& A)
{
    vec x(A.col_dim(), 0.0);
    mul(x, v, A);
    return x;
}

// x += v * a
void acc_mul(vec& x, const vec& v, real_t a);
// x = v * a
void mul(vec& x, const vec& v, real_t a);
inline vec mul(const vec& v, real_t a)
{
    vec x(v.dim(), 0.0);
    mul(x, v, a);
    return x;
}
/* A = Aa;
 */
void scale(matrix& A, real_t a);
void scale(vec& v, real_t a);

sqr_mat part_mat(const matrix& mat,
                 const std::vector<size_t>& sorted_valid_idx);
matrix part_mat(const matrix& mat, const std::vector<size_t>& sorted_row_idx,
                const std::vector<size_t>& sorted_col_idx);
vec part_vec(const vec& v, const std::vector<size_t>& sorted_valid_idx);

sqr_mat transpose(const sqr_mat& mat);
}
