#include "spmatrix_oper.hpp"
#include "spblas.hpp"
namespace markovgg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         spmatrix_format format)
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
    return M.create(format);
}

void dot(vector_mutable_view y, const spmatrix& A, bool transposeA,
         vector_const_view x)
{
    // y = alpha * A * x + beta * y
    spblas_matrix_vector(1.0, A, transposeA, x, 0.0, y);
}

vector dot(const spmatrix& A, bool transposeA, vector_const_view x)
{
    vector y(x.dim());
    dot(y, A, transposeA, x);
    return y;
}
}
