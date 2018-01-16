#pragma once

#include "matrix.hpp"
#include "spmatrix.hpp"
#include "vector.hpp"

namespace scigg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         bool is_row_compressed);

// y = A * x
void dot(spmatrix_const_view A, vector_const_view x, vector_mutable_view y);
inline vector operator*(spmatrix_const_view A, vector_const_view x)
{
    vector y(x.dim());
    dot(A, x, y);
    return y;
}
inline vector operator*(vector_const_view x, spmatrix_const_view A)
{
    vector y(x.dim());
    dot(A.transpose(), x, y);
    return y;
}
inline vector dot(spmatrix_const_view A, vector_const_view x)
{
    vector y(x.dim());
    dot(A, x, y);
    return y;
}

real_t dot(spvec_const_view x, vector_const_view y);
real_t dot(spvec_const_view x, spvec_const_view y);
inline real_t dot(vector_const_view x, spvec_const_view y) { return dot(y, x); }
inline real_t operator*(spvec_const_view x, vector_const_view y)
{
    return dot(x, y);
}
inline real_t operator*(vector_const_view x, spvec_const_view y)
{
    return dot(y, x);
}
inline real_t operator*(spvec_const_view x, spvec_const_view y)
{
    return dot(x, y);
}

// C = A + B
spmatrix add(spmatrix_const_view A, spmatrix_const_view B);
inline spmatrix operator+(spmatrix_const_view A, spmatrix_const_view B)
{
    return add(A, B);
}

matrix spmatrix2dense(spmatrix_const_view A);
spmatrix dense2spmatrix(matrix_const_view A, real_t drop_tol);
}
