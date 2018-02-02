#include <vector>
#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "petrinet.hpp"

TEST(test_petri_net, create)
{
    petri_net_tmpl<double, std::vector<uint_t>, uint_t> pn(2);
    auto t1 =
        pn.add_trans(1, true, 3.0).add_in_arc(0, 1).add_out_arc(1, 1).idx();
    auto t2 =
        pn.add_trans(0, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx();
    auto t3 =
        pn.add_trans(2, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx();
    pn.add_trans(2, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1);
    pn.add_trans(2, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1);
    ASSERT_EQ(pn[t1].prio(), 1);
    ASSERT_EQ(pn[t2].prio(), 0);
    ASSERT_EQ(pn[t3].prio(), 2);
}

TEST(test_petri_net, fire1)
{
    petri_net_tmpl<double, std::vector<uint_t>, uint_t> pn(2);
    auto t1 =
        pn.add_trans(0, true, 3.0).add_in_arc(0, 1).add_out_arc(1, 1).idx();
    auto t2 =
        pn.add_trans(0, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx();
    auto t3 =
        pn.add_trans(0, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx();
    auto init_mk = pn.empty_marking();
    init_mk[0] = 1;
    ASSERT_TRUE(pn.is_trans_enabled(t1, init_mk));
    ASSERT_FALSE(pn.is_trans_enabled(t2, init_mk));
    ASSERT_FALSE(pn.is_trans_enabled(t3, init_mk));
    auto nmk = pn.fire_trans(t1, init_mk);
    ASSERT_EQ(nmk[0], 0);
    ASSERT_EQ(nmk[1], 1);
    ASSERT_FALSE(pn.is_trans_enabled(t1, nmk));
    ASSERT_TRUE(pn.is_trans_enabled(t2, nmk));
    ASSERT_TRUE(pn.is_trans_enabled(t3, nmk));
    auto etrans = pn.enabled_trans(nmk);
    ASSERT_EQ(etrans.size(), 2);
    auto nnmk = pn.fire_trans(etrans[0], nmk);
    ASSERT_EQ(nnmk[0], 1);
    ASSERT_EQ(nnmk[1], 0);
    nnmk = pn.fire_trans(etrans[1], nmk);
    ASSERT_EQ(nnmk[0], 1);
    ASSERT_EQ(nnmk[1], 0);
}

struct mk_comp
{
    bool operator()(const std::vector<uint_t>& v1,
                    const std::vector<uint_t>& v2) const
    {
        for (size_t i = 0; i < v1.size(); i++)
        {
            if (v1[i] < v2[i])
            {
                return true;
            }
            else if (v1[i] > v2[i])
            {
                return false;
            }
        }
        return false;
    }
};

TEST(test_petri_net, gen_reach_graph1)
{
    using petri_net_t = petri_net_tmpl<double, std::vector<uint_t>, uint_t>;
    petri_net_t pn(2);
    auto t1 =
        pn.add_trans(0, true, 1.0).add_in_arc(0, 1).add_out_arc(1, 1).idx();
    auto t2 =
        pn.add_trans(0, true, 2.0).add_in_arc(1, 1).add_out_arc(0, 1).idx();
    auto init_mk = pn.empty_marking();
    init_mk[0] = 1;
    using reach_graph = reach_graph_tmpl<petri_net_t, ordered_digraph, mk_comp>;
    auto rg = reach_graph();
    gen_reach_graph(pn, init_mk, rg);
    ASSERT_EQ(rg.node_count(), 2);
    ASSERT_EQ(rg[0].data()[0], 1);
    ASSERT_EQ(rg[0].data()[1], 0);
    ASSERT_EQ(rg[0].out_arc().size(), 1);
    ASSERT_EQ(rg[0].out_arc()[0].data()->idx(), t1);
    ASSERT_EQ(rg[0].out_arc()[0].data()->data(), 1.0);
    ASSERT_EQ(rg[1].data()[0], 0);
    ASSERT_EQ(rg[1].data()[1], 1);
    ASSERT_EQ(rg[1].out_arc().size(), 1);
    ASSERT_EQ(rg[1].out_arc()[0].data()->idx(), t2);
    ASSERT_EQ(rg[1].out_arc()[0].data()->data(), 2.0);
}

TEST(test_petri_net, gen_reach_graph2)
{
    using petri_net_t = default_petri_net<double>;
    petri_net_t pn(2);
    auto t1 =
        pn.add_trans(0, true, 1.0).add_in_arc(0, 1).add_out_arc(1, 1).idx();
    auto t2 =
        pn.add_trans(0, true, 2.0).add_in_arc(1, 1).add_out_arc(0, 1).idx();
    auto init_mk = pn.empty_marking();
    init_mk[0] = 1;
    using reach_graph = default_reach_graph<petri_net_t>;
    auto rg = reach_graph();
    gen_reach_graph(pn, init_mk, rg);
    ASSERT_EQ(rg.node_count(), 2);
    ASSERT_EQ(rg[0].data()[0], 1);
    ASSERT_EQ(rg[0].data()[1], 0);
    ASSERT_EQ(rg[0].out_arc().size(), 1);
    ASSERT_EQ(rg[0].out_arc()[0].data()->idx(), t1);
    ASSERT_EQ(rg[0].out_arc()[0].data()->data(), 1.0);
    ASSERT_EQ(rg[1].data()[0], 0);
    ASSERT_EQ(rg[1].data()[1], 1);
    ASSERT_EQ(rg[1].out_arc().size(), 1);
    ASSERT_EQ(rg[1].out_arc()[0].data()->idx(), t2);
    ASSERT_EQ(rg[1].out_arc()[0].data()->data(), 2.0);
}
