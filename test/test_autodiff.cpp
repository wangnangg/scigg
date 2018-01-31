#include <memory>
#include "autodiff.hpp"
#include "common.hpp"
#include "gtest/gtest.h"

TEST(test_autodiff, diff1_sum)
{
    // s = x + y + z
    auto st = diff1_expr_store();
    auto s = st.sum(st.var(0), st.var(1), st.var(2));

    const auto val = create_vector(3, {1.0, 2.0, 3.0});
    auto grad = vector(3);
    auto result = s->eval(val, grad);
    ASSERT_EQ(result, 6.0);
    ASSERT_EQ(grad, create_vector(3, {1.0, 1.0, 1.0}));
}

TEST(test_autodiff, diff1_prod)
{
    // s = (x + y ) * (x + z) * 5
    auto st = diff1_expr_store();
    auto x = st.var(0);
    auto y = st.var(1);
    auto z = st.var(2);
    auto s = st.prod(st.sum(x, y), st.sum(x, z), st.con(5));

    const auto val = create_vector(3, {1, 2, 3});
    auto grad = vector(3);
    auto result = s->eval(val, grad);
    ASSERT_EQ(result, 60.0);
    ASSERT_EQ(grad, create_vector(3, {35, 20, 15}));
}
