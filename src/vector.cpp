#include "vector.hpp"

namespace markovgg
{
bool operator==(const vector &x, const vector &y)
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
