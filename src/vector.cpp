#include "vector.hpp"

namespace markovgg
{
vector_mutable_view::vector_mutable_view(vector &v)
    : vector_mutable_view(&v[0], v.dim(), v.inc())
{
}

vector_const_view::vector_const_view(const vector &v)
    : vector_const_view(&v[0], v.dim(), v.inc())
{
}

vector_const_view::vector_const_view(vector_mutable_view v)
    : vector_const_view(&v[0], v.dim(), v.inc())
{
}

bool operator==(vector_const_view x, vector_const_view y)
{
    if (x.dim() != y.dim())
    {
        return false;
    }
    for (size_t i = 0; i < x.dim(); i++)
    {
        if (x[i] != y[i])
        {
            return false;
        }
    }
    return true;
}
vector create_vector(size_t n, const std::vector<double> &v)
{
    vector vo(n);
    for (size_t i = 0; i < n; i++)
    {
        vo[i] = v[i];
    }
    return vo;
}
}
