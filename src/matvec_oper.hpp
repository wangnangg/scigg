#pragma once
#include "blas.hpp"
#include "matrix.hpp"
#include "vector.hpp"

namespace markovgg
{
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
void fill(matrix_mutable_view A, real_t val);

// y = x
void copy(vector_mutable_view y, vector_const_view x);
void copy(matrix_mutable_view Y, matrix_const_view X);

inline void scale(vector_mutable_view x, real_t alpha) { blas_scale(alpha, x); }
inline vector_mutable_view operator*=(vector_mutable_view x, real_t alpha)
{
    scale(x, alpha);
    return x;
}

vector operator*(real_t alpha, vector_const_view x);
inline vector operator*(vector_const_view x, real_t alpha) { return alpha * x; }

void scale(matrix_mutable_view A, real_t alpha);
inline matrix_mutable_view operator*=(matrix_mutable_view A, real_t alpha)
{
    scale(A, alpha);
    return A;
}

matrix operator*(real_t alpha, matrix_const_view A);
inline matrix operator*(matrix_const_view A, real_t alpha) { return alpha * A; }

inline real_t norm1(vector_const_view vec) { return blas_abs_sum(vec); }

inline real_t norm2(vector_const_view vec) { return blas_norm2(vec); }

inline void set_norm1(vector_mutable_view vec, real_t val)
{
    real_t n1 = norm1(vec);
    scale(vec, val / n1);
}
inline void set_norm2(vector_mutable_view vec, real_t val)
{
    real_t n2 = norm2(vec);
    scale(vec, val / n2);
}

// add z = x + y
void add(vector_mutable_view z, vector_const_view x, vector_const_view y);
inline vector add(vector_const_view x, vector_const_view y)
{
    vector z(x.dim());
    add(z, x, y);
    return z;
}
inline vector operator+(vector_const_view x, vector_const_view y)
{
    return add(x, y);
}

// inc z += x
void inc(vector_mutable_view z, vector_const_view x);
inline vector_mutable_view operator+=(vector_mutable_view z,
                                      vector_const_view x)
{
    inc(z, x);
    return z;
}

// add z = x - y
void sub(vector_mutable_view z, vector_const_view x, vector_const_view y);
inline vector sub(vector_const_view x, vector_const_view y)
{
    vector z(x.dim());
    sub(z, x, y);
    return z;
}
inline vector operator-(vector_const_view x, vector_const_view y)
{
    return sub(x, y);
}

// z -= x
void dec(vector_mutable_view z, vector_const_view x);
inline vector_mutable_view operator-=(vector_mutable_view z,
                                      vector_const_view x)
{
    dec(z, x);
    return z;
}

// dot x * y
real_t dot(vector_const_view x, vector_const_view y);
inline real_t operator*(vector_const_view x, vector_const_view y)
{
    return dot(x, y);
}

// dot z = A . x
void dot(vector_mutable_view z, matrix_const_view A, vector_const_view x);
inline vector dot(matrix_const_view A, vector_const_view x)
{
    vector z(A.m());
    dot(z, A, x);
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
void dot(matrix_mutable_view C, matrix_const_view A, matrix_const_view B);
inline matrix dot(matrix_const_view A, matrix_const_view B)
{
    matrix C(A.m(), B.n());
    dot(C, A, B);
    return C;
}
inline matrix operator*(matrix_const_view A, matrix_const_view B)
{
    return dot(A, B);
}
}
