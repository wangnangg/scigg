#pragma once
#include <vector>
#include "spmatvec.hpp"

namespace scigg
{
struct scc_component
{
    std::vector<size_t> node_list;
    bool is_bottom;
};

// digraph must be in from_to format.
std::vector<scc_component> strongly_connected(spmatrix_const_view digraph);
}
