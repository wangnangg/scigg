#pragma once

#include "spmatrix.hpp"

namespace markovgg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         spmatrix_format format = CPR_ROW);
}
