#pragma once
#include <vector>
#include "spmatrix.hpp"

namespace markovgg
{
struct scc_component
{
    std::vector<size_t> node_list;
    bool is_bottom;
};

// digraph must be in from_to format.
std::vector<scc_component> strongly_connected(const spmatrix& digraph);
}
