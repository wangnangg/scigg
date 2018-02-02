#pragma once
#include <set>
#include <vector>
#include "type.hpp"
namespace scigg
{
// data can be compared
template <typename Tn, typename Ta, typename TnCompare = std::less<Tn>>
class ordered_digraph
{
public:
    using nid_t = size_t;
    class arc
    {
        nid_t _dst;
        Ta _data;

    public:
        arc(nid_t dst, Ta data) : _dst(dst), _data(std::move(data)) {}
        const Ta& data() const { return _data; }
    };
    class node
    {
        nid_t _idx;
        Tn _data;
        std::vector<arc> _out_arc;

    public:
        node(nid_t idx, Tn data) : _idx(idx), _data(std::move(data)) {}
        node& add_arc(nid_t dst, Ta data)
        {
            _out_arc.emplace_back(dst, std::move(data));
            return *this;
        }
        nid_t idx() const { return _idx; }
        const Tn& data() const { return _data; }
        const std::vector<arc>& out_arc() const { return _out_arc; }
    };
    using ndata_comp = bool(const Tn& d1, const Tn& d2);

private:
    struct node_comp
    {
        const TnCompare& ndata_comp;
        node_comp(const TnCompare& comp) : ndata_comp(comp) {}
        bool operator()(const node& n1, const node& n2) const
        {
            return ndata_comp(n1.data(), n2.data());
        }
    };
    std::set<node, node_comp> _node_set;
    std::vector<node*> _nid_lookup;

public:  // observer
    const node& operator[](nid_t nid) const { return *_nid_lookup[nid]; }
    node& operator[](nid_t nid) { return *_nid_lookup[nid]; }
    nid_t node_count() const { return _node_set.size(); }

public:  // modifier
    ordered_digraph(const TnCompare& comp = TnCompare())
        : _node_set(node_comp(comp))
    {
    }
    node& add_node(Tn data)
    {
        nid_t new_idx = node_count();
        auto pair = _node_set.emplace(new_idx, std::move(data));
        node* node_ptr = const_cast<node*>(&*(pair.first));
        if (pair.second)
        {  // new node inserted
            _nid_lookup.push_back(node_ptr);
        }
        return *node_ptr;
    }
    node& add_arc(nid_t src, nid_t dst, Ta data)
    {
        assert(src < node_count());
        return (*this)[src].add_arc(dst, std::move(data));
    }
};
}
