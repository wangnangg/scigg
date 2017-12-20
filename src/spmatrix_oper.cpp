#include "spmatrix_oper.hpp"
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
}
