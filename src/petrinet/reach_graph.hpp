#pragma once
#include "graph.hpp"
#include "petri_net.hpp"
#include "type.hpp"
namespace scigg
{
template <typename PetriNet,
          template <typename Tn, typename Ta, typename TnFunc> typename DiGraph,
          typename MkFunc>
using reach_graph_tmpl = DiGraph<typename PetriNet::marking,
                                 const typename PetriNet::transition*, MkFunc>;

template <typename PetriNet,
          template <typename Tn, typename Ta, typename TnFunc> typename DiGraph,
          typename MkFunc>
void gen_reach_graph(const PetriNet& pn,
                     const typename PetriNet::marking& init_mk,
                     reach_graph_tmpl<PetriNet, DiGraph, MkFunc>& rg)
{
    rg.add_node(init_mk);
    typename reach_graph_tmpl<PetriNet, DiGraph, MkFunc>::nid_t curr_nid = 0;
    while (curr_nid < rg.node_count())
    {
        const typename PetriNet::marking& curr_mk = rg[curr_nid].data();
        for (const auto& tid : pn.enabled_trans(curr_mk))
        {
            auto new_mk = pn.fire_trans(tid, curr_mk);
            size_t dst_nid = rg.add_node(std::move(new_mk)).idx();
            rg.add_arc(curr_nid, dst_nid, &pn[tid]);
        }
        curr_nid += 1;
    }
}
}
