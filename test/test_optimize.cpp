#include <cmath>
#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "linalg.hpp"
#include "matvec.hpp"
#include "optimize.hpp"

// y = x0 * x0 + 2 * (x1 - 1) * (x1 - 1) + 5
// dy/dx0 = 2 * x0
// dy/dx1 = 4 * (x1 -  1)
// y* = 5, x* = (0, 1);

real_t obj_function1(vector_const_view x, vector_mutable_view grad)
{
    std::cout << "trying: " << x << std::endl;
    real_t y = x[0] * x[0] + 2.0 * (x[1] - 1) * (x[1] - 1) + 5.0;
    grad[0] = 2 * x[0];
    grad[1] = 4 * (x[1] - 1);
    return y;
}

real_t obj_wood(vector_const_view x, vector_mutable_view grad)
{
    std::cout << "trying: " << x << std::endl;
    real_t y = 100 * std::pow(x[0] * x[0] - x[1], 2) +                     //
               std::pow(x[0] - 1, 2) +                                     //
               std::pow(x[2] - 1, 2) +                                     //
               90 * std::pow(x[2] * x[2] - x[3], 2) +                      //
               10.1 * (std::pow(x[1] - 1, 2) + std::pow(x[3] - 1.0, 2)) +  //
               19.8 * (x[1] - 1) * (x[3] - 1);
    grad[0] = 200 * (x[0] * x[0] - x[1]) * 2 * x[0] + 2 * (x[0] - 1);
    grad[1] =
        -200 * (x[0] * x[0] - x[1]) + 20.2 * (x[1] - 1) + 19.8 * (x[3] - 1);
    grad[2] = 2 * (x[2] - 1) + 180 * (x[2] * x[2] - x[3]) * 2 * x[2];
    grad[3] =
        -180 * (x[2] * x[2] - x[3]) + 20.2 * (x[3] - 1.0) + 19.8 * (x[1] - 1);
    return y;
}

TEST(test_optimize, quasi_newton_bfgs1)
{
    auto x = create_vector(2, {10, 10});
    real_t y;
    real_t tol = 1e-10;
    uint_t max_iter = 100;
    real_t gradn =
        quasi_newton_bfgs(obj_function1, x, y, 0.1, 1e-4, 0.9, tol, max_iter);
    ASSERT_LT(gradn, tol);
    ASSERT_NEAR(y, 5.0, 1e-6);
    ASSERT_TRUE(near_eq(x, create_vector(2, {0, 1}), 1e-6));
}

TEST(test_optimize, quasi_newton_bfgs2)
{
    auto x = create_vector(4, {-3, -1, -3, -1});
    real_t y;
    real_t tol = 1e-10;
    uint_t max_iter = 100;
    real_t gradn =
        quasi_newton_bfgs(obj_wood, x, y, 0.1, 1e-4, 0.9, tol, max_iter);
    ASSERT_LT(gradn, tol);
    ASSERT_NEAR(y, 0.0, 1e-6);
    ASSERT_TRUE(near_eq(x, create_vector(4, {1, 1, 1, 1}), 1e-6));
}

TEST(test_optimize, obj_wood)
{
    auto x = create_vector(4, {1, 1, 1, 1});
    vector grad(x.dim());
    real_t y;
    y = obj_wood(x, grad);
    ASSERT_DOUBLE_EQ(y, 0.0);
    ASSERT_TRUE(near_eq(grad, create_vector(4, {0, 0, 0, 0}), 1e-10));
    x = create_vector(4, {1, 2, 3, 4});
    y = obj_wood(x, grad);
    ASSERT_DOUBLE_EQ(y, 2514.4);
    print(grad);
    ASSERT_TRUE(
        near_eq(grad, create_vector(4, {-400, 279.6, 5404, -819.6}), 1e-10));
}
