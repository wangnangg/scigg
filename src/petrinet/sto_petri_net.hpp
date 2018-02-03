#pragma once
#include <functional>
#include "autodiff.hpp"
#include "markov.hpp"
#include "petri_net.hpp"
#include "reach_graph.hpp"

namespace scigg
{
template <typename Marking>
using diff1_spn = petri_net_tmpl<const diff1_expr*, Marking>;

using default_diff1_spn = diff1_spn<default_marking>;

using default_diff1_reach_graph = default_reach_graph<default_diff1_spn>;

template <typename ReachGraph>
std::vector<spmatrix> diff1_reach2markov(const ReachGraph& rg,
                                         vector_const_view param);

template <typename ReachGraph>
class diff1_spn_ss_prob_reward : public diff1_expr
{
public:
    class context
    {
    private:
        const ReachGraph& _rg;
        const typename ReachGraph::node& _node;

    public:
        context(const ReachGraph& rg, size_t nid) : _rg(rg), _node(rg[nid]) {}
        auto ntoken(pid_t pid) { return _node.data()[pid]; }
    };
    using reward_func =
        std::function<real_t(context ct, vector_mutable_view grad)>;

private:
    const ReachGraph& _rg;
    vector_const_view _init_prob;
    matrix_const_view _init_prob_jaco;
    reward_func _rfunc;

public:
    diff1_spn_ss_prob_reward(const ReachGraph& rg, vector_const_view init_prob,
                             matrix_const_view init_prob_jaco,
                             reward_func rfunc)
        : _rg(rg),
          _init_prob(init_prob),
          _init_prob_jaco(init_prob_jaco),
          _rfunc(rfunc)
    {
    }
    real_t eval(vector_const_view x, vector_mutable_view grad) const override
    {
        auto mat_list = diff1_reach2markov(_rg, x);
        const spmatrix& Q = mat_list[0];
        size_t N = _rg.node_count();
        vector ss_prob(N);
        general_markov_ss_prob(Q, _init_prob, ss_prob);
        matrix ss_prob_grad(N, x.dim());
        for (size_t i = 0; i < x.dim(); i++)
        {
            // for param x[i]
            const spmatrix& dQ = mat_list[i + 1];
            general_markov_ss_dprob(Q, _init_prob, dQ, _init_prob_jaco.col(i),
                                    ss_prob, ss_prob_grad.col(i));
        }
        return 0.0;
    }
};
}
