#include <vector>
#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "petrinet.hpp"

TEST(test_petri_net, create)
{
    petri_net_tmpl<double, default_marking> pn(2);
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
    petri_net_tmpl<double, default_marking> pn(2);
    auto t1 = pn.add_trans(0, true, [](auto) { return 3.0; })
                  .add_in_arc(0, 1)
                  .add_out_arc(1, 1)
                  .idx();
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

TEST(test_petri_net, gen_reach_graph1)
{
    using petri_net_t = petri_net_tmpl<double, default_marking>;
    petri_net_t pn(2);
    pn.add_trans(0, true, 1.0).add_in_arc(0, 1).add_out_arc(1, 1).idx();
    pn.add_trans(0, true, 2.0).add_in_arc(1, 1).add_out_arc(0, 1).idx();
    auto init_mk = pn.empty_marking();
    init_mk[0] = 1;
    using reach_graph =
        reach_graph_tmpl<petri_net_t, ordered_digraph, default_marking_comp>;
    auto rg = reach_graph();
    gen_reach_graph(pn, init_mk, rg);
    ASSERT_EQ(rg.node_count(), 2);
    ASSERT_EQ(rg[0].data()[0], 1);
    ASSERT_EQ(rg[0].data()[1], 0);
    ASSERT_EQ(rg[0].out_arc().size(), 1);
    ASSERT_EQ(rg[0].out_arc()[0].data(), 1.0);
    ASSERT_EQ(rg[1].data()[0], 0);
    ASSERT_EQ(rg[1].data()[1], 1);
    ASSERT_EQ(rg[1].out_arc().size(), 1);
    ASSERT_EQ(rg[1].out_arc()[0].data(), 2.0);
}
TEST(test_petri_net, gen_reach_graph2)
{
    using petri_net_t = default_petri_net<double>;
    petri_net_t pn(2);
    pn.add_trans(0, true, [=](auto) { return 1.0; })
        .add_in_arc(0, 1)
        .add_out_arc(1, 1)
        .idx();
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
    ASSERT_EQ(rg[0].out_arc()[0].data(), 1.0);
    ASSERT_EQ(rg[1].data()[0], 0);
    ASSERT_EQ(rg[1].data()[1], 1);
    ASSERT_EQ(rg[1].out_arc().size(), 1);
    ASSERT_EQ(rg[1].out_arc()[0].data(), 2.0);
}

TEST(test_petri_net, diff1_spn)
{
    default_diff1_spn spn(2);
    diff1_expr_store st;
    auto t1 = spn.add_trans(0, true,
                            [&](auto ct) {
                                return st.prod(st.con(ct.ntoken(0)), st.var(0));
                            })
                  .add_in_arc(0)
                  .add_out_arc(1)
                  .idx();
    auto t2 = spn.add_trans(0, true,
                            [&](auto ct) {
                                return st.prod(st.con(ct.ntoken(1)), st.var(0));
                            })
                  .add_in_arc(1, 1)
                  .add_out_arc(0, 1)
                  .idx();
    auto init_mk = spn.empty_marking();
    init_mk[0] = 2;
    auto parm = create_vector(1, {3});
    auto grad = vector(1, 0.0);
    ASSERT_EQ(spn.trans_data(t1, init_mk)->eval(parm, grad), 6.0);
    ASSERT_EQ(grad, create_vector(1, {2.0}));
    ASSERT_EQ(spn.trans_data(t2, init_mk)->eval(parm, grad), 0.0);
    ASSERT_EQ(grad, create_vector(1, {0.0}));
}
