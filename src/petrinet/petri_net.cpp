#include "petri_net.hpp"

namespace scigg
{
bool default_marking_comp::operator()(const default_marking& m1,
                                      const default_marking& m2) const
{
    for (size_t i = 0; i < m1.size(); i++)
    {
        if (m1[i] < m2[i])
        {
            return true;
        }
        else if (m1[i] > m2[i])
        {
            return false;
        }
    }
    return false;
}
}
