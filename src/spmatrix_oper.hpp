#pragma once

#include "matrix.hpp"
#include "spmatrix.hpp"
#include "vector.hpp"

namespace markovgg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         bool is_row_compressed);

// y = A * x
void dot(vector_mutable_view y, spmatrix_const_view A, vector_const_view x);
inline vector operator*(spmatrix_const_view A, vector_const_view x)
{
    vector y(x.dim());
    dot(y, A, x);
    return y;
}
inline vector operator*(vector_const_view x, spmatrix_const_view A)
{
    vector y(x.dim());
    dot(y, A.transpose(), x);
    return y;
}
inline vector dot(spmatrix_const_view A, vector_const_view x)
{
    vector y(x.dim());
    dot(y, A, x);
    return y;
}

matrix spmatrix2dense(spmatrix_const_view A);
}
