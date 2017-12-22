#include "digraph_algo.hpp"
namespace markovgg
{
struct scc_node
{
    size_t depth;
    size_t lowlink;
    bool onstack;
    scc_node() : depth(0), lowlink(0), onstack(false) {}
};

bool recursive_visit(const spmatrix& mat, std::vector<scc_node>& node_info,
                     size_t src_idx, size_t& depth,
                     std::vector<size_t>& node_stack,
                     std::vector<scc_component>& scc_list)
{
    auto& src_node = node_info[src_idx];
    src_node.depth = depth;
    src_node.lowlink = depth;
    src_node.onstack = true;
    depth += 1;
    node_stack.push_back(src_idx);
    bool is_bottom = true;
    for (const auto& e : mat[src_idx])
    {
        size_t dst_idx = e.idx;
        if (dst_idx != src_idx)
        {
            auto& dst_node_info = node_info[dst_idx];
            if (dst_node_info.depth == 0)
            {
                if (!recursive_visit(mat, node_info, dst_idx, depth, node_stack,
                                     scc_list))
                {
                    is_bottom = false;
                }
                if (src_node.lowlink > dst_node_info.lowlink)
                {
                    src_node.lowlink = dst_node_info.lowlink;
                }
            }
            else if (dst_node_info.onstack)
            {
                if (src_node.lowlink > dst_node_info.depth)
                {
                    src_node.lowlink = dst_node_info.depth;
                }
            }
            else
            {
                is_bottom = false;
            }
        }
    }

    if (src_node.lowlink == src_node.depth)
    {
        std::vector<size_t> node_list;
        size_t top_node_idx;
        do
        {
            top_node_idx = node_stack.back();
            node_stack.pop_back();
            node_list.push_back(top_node_idx);
            node_info[top_node_idx].onstack = false;
        } while (top_node_idx != src_idx);
        scc_list.push_back({std::move(node_list), is_bottom});
        return false;
    }
    else
    {
        return is_bottom;
    }
}

std::vector<scc_component> strongly_connected(const spmatrix& digraph)
{
    assert(digraph.m() == digraph.n());
    std::vector<scc_component> scc_list;
    std::vector<scc_node> node_info(digraph.m());
    std::vector<size_t> scc_stack;
    size_t depth = 1;
    for (size_t i = 0; i < digraph.m(); i++)
    {
        if (node_info[i].depth == 0)
        {
            recursive_visit(digraph, node_info, i, depth, scc_stack, scc_list);
        }
    }
    return scc_list;
}
}
