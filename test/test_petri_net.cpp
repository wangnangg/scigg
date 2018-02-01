#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "petrinet.hpp"

TEST(test_petri_net, create)
{
    petri_net<double, std::vector<uint_t>, uint_t> pn(2);
    auto t1 = pn.add_trans(1, true, 3.0).add_in_arc(0, 1).add_out_arc(1, 1).idx;
    auto t2 = pn.add_trans(0, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx;
    auto t3 = pn.add_trans(2, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx;
    pn.add_trans(2, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1);
    pn.add_trans(2, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1);
    ASSERT_EQ(pn[t1].prio, 1);
    ASSERT_EQ(pn[t2].prio, 0);
    ASSERT_EQ(pn[t3].prio, 2);
}

TEST(test_petri_net, fire1)
{
    petri_net<double, std::vector<uint_t>, uint_t> pn(2);
    auto t1 = pn.add_trans(0, true, 3.0).add_in_arc(0, 1).add_out_arc(1, 1).idx;
    auto t2 = pn.add_trans(0, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx;
    auto t3 = pn.add_trans(0, true, 3.0).add_in_arc(1, 1).add_out_arc(0, 1).idx;
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
