#include "spmatrix_oper.hpp"
#include "spblas.hpp"
namespace markovgg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         bool is_row_compressed)
{
    assert(m * n == v.size());
    spmatrix_creator M(m, n);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            M.add_entry(i, j, v[i * n + j]);
        }
    }
    return M.create(is_row_compressed);
}

void dot(vector_mutable_view y, spmatrix_const_view A, vector_const_view x)
{
    // y = alpha * A * x + beta * y
    spblas_matrix_vector(1.0, A, x, 0.0, y);
}

vector dot(spmatrix_const_view A, vector_const_view x)
{
    vector y(x.dim());
    dot(y, A, x);
    return y;
}
}
