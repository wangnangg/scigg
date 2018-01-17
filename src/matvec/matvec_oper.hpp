#pragma once
#include "blas.hpp"
#include "matrix.hpp"
#include "vector.hpp"

namespace scigg
{
inline real_t abs(real_t v) { return v >= 0 ? v : -v; }
inline real_t max(real_t x, real_t y) { return x > y ? x : y; }
inline real_t min(real_t x, real_t y) { return x > y ? y : x; }

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

void fill(real_t val, vector_mutable_view x);
void fill(real_t val, matrix_mutable_view A);

// y = x
void copy(vector_const_view x, vector_mutable_view y);
void copy(matrix_const_view X, matrix_mutable_view Y);

inline void scale(real_t alpha, vector_mutable_view x) { blas_scale(alpha, x); }
inline vector_mutable_view operator*=(vector_mutable_view x, real_t alpha)
{
    scale(alpha, x);
    return x;
}

vector operator*(real_t alpha, vector_const_view x);
inline vector operator*(vector_const_view x, real_t alpha) { return alpha * x; }

void scale(real_t alpha, matrix_mutable_view A);
inline matrix_mutable_view operator*=(matrix_mutable_view A, real_t alpha)
{
    scale(alpha, A);
    return A;
}

matrix operator*(real_t alpha, matrix_const_view A);
inline matrix operator*(matrix_const_view A, real_t alpha) { return alpha * A; }

inline real_t norm1(vector_const_view vec) { return blas_abs_sum(vec); }

inline real_t norm2(vector_const_view vec) { return blas_norm2(vec); }

inline void set_norm1(real_t val, vector_mutable_view vec)
{
    real_t n1 = norm1(vec);
    scale(val / n1, vec);
}
inline void set_norm2(real_t val, vector_mutable_view vec)
{
    real_t n2 = norm2(vec);
    scale(val / n2, vec);
}

// add z = x + y
void add(vector_const_view x, vector_const_view y, vector_mutable_view z);
inline vector add(vector_const_view x, vector_const_view y)
{
    vector z(x.dim());
    add(x, y, z);
    return z;
}
inline vector operator+(vector_const_view x, vector_const_view y)
{
    return add(x, y);
}

void add(matrix_const_view x, matrix_const_view y, matrix_mutable_view z);
inline matrix add(matrix_const_view x, matrix_const_view y)
{
    auto z = matrix(x.m(), x.n());
    add(x, y, z);
    return z;
}
inline matrix operator+(matrix_const_view x, matrix_const_view y)
{
    return add(x, y);
}

// inc z += x
void inc(vector_const_view x, vector_mutable_view z);
inline vector_mutable_view operator+=(vector_mutable_view z,
                                      vector_const_view x)
{
    inc(x, z);
    return z;
}

// add z = x - y
void sub(vector_const_view x, vector_const_view y, vector_mutable_view z);
inline vector sub(vector_const_view x, vector_const_view y)
{
    vector z(x.dim());
    sub(x, y, z);
    return z;
}
inline vector operator-(vector_const_view x, vector_const_view y)
{
    return sub(x, y);
}

// z -= x
void dec(vector_const_view x, vector_mutable_view z);
inline vector_mutable_view operator-=(vector_mutable_view z,
                                      vector_const_view x)
{
    dec(x, z);
    return z;
}

// dot x * y
real_t dot(vector_const_view x, vector_const_view y);
inline real_t operator*(vector_const_view x, vector_const_view y)
{
    return dot(x, y);
}

// dot z = A . x
void dot(matrix_const_view A, vector_const_view x, vector_mutable_view z);
inline vector dot(matrix_const_view A, vector_const_view x)
{
    vector z(A.m());
    dot(A, x, z);
    return z;
}
inline vector operator*(matrix_const_view A, vector_const_view x)
{
    return dot(A, x);
}
inline vector operator*(vector_const_view x, matrix_const_view A)
{
    return dot(A.transpose(), x);
}

// dot C = A . B
void dot(matrix_const_view A, matrix_const_view B, matrix_mutable_view C);
inline matrix dot(matrix_const_view A, matrix_const_view B)
{
    matrix C(A.m(), B.n());
    dot(A, B, C);
    return C;
}
inline matrix operator*(matrix_const_view A, matrix_const_view B)
{
    return dot(A, B);
}
}
