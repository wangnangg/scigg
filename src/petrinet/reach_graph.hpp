#pragma once
#include <memory>
#include "graph.hpp"
#include "petri_net.hpp"
#include "type.hpp"

namespace scigg
{
template <typename PetriNet,
          template <typename Tn, typename Ta, typename TnFunc> typename DiGraph,
          typename MkFunc>
using reach_graph_tmpl =
    DiGraph<typename PetriNet::marking, typename PetriNet::tdata, MkFunc>;

template <typename PetriNet>
using default_reach_graph =
    reach_graph_tmpl<PetriNet, ordered_digraph, default_marking_comp>;

template <typename PetriNet, typename ReachGraph>
void gen_reach_graph(const PetriNet& pn,
                     const typename PetriNet::marking& init_mk, ReachGraph& rg)
{
    rg.add_node(init_mk);
    typename ReachGraph::nid_t curr_nid = 0;
    while (curr_nid < rg.node_count())
    {
        const typename PetriNet::marking& curr_mk = rg[curr_nid].data();
        for (const auto& tid : pn.enabled_trans(curr_mk))
        {
            auto new_mk = pn.fire_trans(tid, curr_mk);
            size_t dst_nid = rg.add_node(std::move(new_mk)).idx();
            rg.add_arc(curr_nid, dst_nid, pn.trans_data(tid, curr_mk));
        }
        curr_nid += 1;
    }
}
}
