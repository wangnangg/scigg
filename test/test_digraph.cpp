#include "common.hpp"
#include "debug_utils.hpp"
#include "digraph_algo.hpp"
#include "gtest/gtest.h"

size_t bottom_count(const std::vector<scc_component>& scc_list)
{
    size_t counter = 0;
    for (const auto& scc : scc_list)
    {
        if (scc.is_bottom)
        {
            counter += 1;
        }
    }
    return counter;
}

TEST(test_digraph, strongly_connected_components)
{
    spmatrix_creator crtor(10, 10);
    crtor.add_entry(0, 1, 1);
    crtor.add_entry(1, 0, 1);
    crtor.add_entry(1, 3, 1);
    crtor.add_entry(3, 6, 1);
    crtor.add_entry(3, 4, 1);
    crtor.add_entry(4, 5, 1);
    crtor.add_entry(4, 9, 1);
    crtor.add_entry(5, 3, 1);
    crtor.add_entry(5, 8, 1);
    crtor.add_entry(6, 7, 1);
    crtor.add_entry(7, 6, 1);
    crtor.add_entry(8, 9, 1);
    crtor.add_entry(9, 8, 1);
    const auto scc_list = strongly_connected(crtor.create(true));
    ASSERT_EQ(bottom_count(scc_list), 3);
}
