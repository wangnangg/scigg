#pragma once
#include <vector>
#include "type.hpp"

template <typename Tn, typename Ta>
class digraph
{
public:
    struct arc
    {
        const size_t dst;
        Ta data;
        arc(size_t dst, Ta data) : dst(dst), data(std::move(data)) {}
    };
    struct node
    {
        const size_t idx;
        Tn data;
        std::vector<arc> out_arc;
        node(size_t idx, Tn data) : idx(idx), data(std::move(data)) {}
        node& add_arc(size_t dst, Ta data)
        {
            out_arc.emplace_back(dst, std::move(data));
            return *this;
        }
    };

private:
    std::vector<node> _node_list;

public:  // observer
    const node& operator[](size_t nid) const { return _node_list[nid]; }
    node& operator[](size_t nid) { return _node_list[nid]; }
    size_t node_count() const { return _node_list.size(); }

public:  // modifier
    node& add_node(Tn data)
    {
        size_t new_idx = _node_list.size();
        _node_list.emplace_back(new_idx, std::move(data));
        return _node_list.back();
    }
    node& add_arc(size_t src, size_t dst, Ta data)
    {
        assert(src < node_count());
        return _node_list[src].add_arc(dst, std::move(data));
    }
};
