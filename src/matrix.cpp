#include "matrix.hpp"

namespace markovgg
{
bool operator==(matrix_const_view m1, matrix_const_view m2)
{
    if (m1.m() != m2.m() || m1.n() != m2.n())
    {
        return false;
    }
    for (size_t i = 0; i < m1.m(); i++)
    {
        for (size_t j = 0; j < m1.n(); j++)
        {
            if (m1(i, j) != m2(i, j))
            {
                return false;
            }
        }
    }
    return true;
}
}
