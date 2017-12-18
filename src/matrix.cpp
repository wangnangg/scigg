#include "matrix.hpp"

namespace markovgg
{
bool operator==(const matrix &m1, const matrix &m2)
{
    if (m1.m() != m2.m() || m1.n() != m2.n())
    {
        return false;
    }
    const real_t *e1 = m1.begin();
    const real_t *e2 = m2.begin();
    while (e1 != m1.end())
    {
        if (*e1 != *e2)
        {
            return false;
        }
        e1 += 1;
        e2 += 1;
    }
    return true;
}
}
