#include "vector.hpp"

namespace markovgg
{
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
}
