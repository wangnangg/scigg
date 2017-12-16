#include <iostream>
#include <vector>
#include "blas.hpp"
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"

TEST(test_blas, level1)
{
    vector v1(3, 0);
    v1[0] = 1;
    v1[1] = 2;
    v1[2] = 3;
    vector v2(3, 0);
    v2[0] = 4;
    v2[1] = 5;
    v2[2] = 6;

    ASSERT_EQ(dot(v1, v2), 32);

    vector v3(2, 0);
    v3[0] = 3;
    v3[1] = 4;
    ASSERT_EQ(norm2(v3), 5);

    ASSERT_EQ(abs_sum(v1), 6);

    ASSERT_EQ(abs_max_idx(v2), 2);

    swap(v1, v2);
    ASSERT_EQ(v1[0], 4);
    ASSERT_EQ(v1[1], 5);
    ASSERT_EQ(v1[2], 6);
    swap(v1, v2);

    vector v4(3, 0);
    copy(v1, v4);
    ASSERT_TRUE(v4 == v1);

    axpy(2.0, v1, v2);
    ASSERT_EQ(v2[0], 6);
    ASSERT_EQ(v2[1], 9);
    ASSERT_EQ(v2[2], 12);
    axpy(-2.0, v1, v2);

    scale(10, v1);
    ASSERT_EQ(v1[0], 10);
    ASSERT_EQ(v1[1], 20);
    ASSERT_EQ(v1[2], 30);
    scale(0.1, v1);
}

matrix create_matrix(size_t m, size_t n, const std::vector<double>& v)
{
    matrix M(m, n);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            M(i, j) = v[i * n + j];
        }
    }
    return M;
}
vector create_vector(size_t n, const std::vector<double>& v)
{
    vector vo(n);
    for (size_t i = 0; i < n; i++)
    {
        vo[i] = v[i];
    }
    return vo;
}
TEST(test_blas, level2)
{
    matrix M = create_matrix(4, 5,
                             {
                                 1,  2,  3,  4,  5,   //
                                 6,  7,  8,  9,  10,  //
                                 11, 12, 13, 14, 15,  //
                                 16, 17, 18, 19, 20   //
                             });

    vector x = create_vector(5, {5, 4, 3, 2, 1});
    vector y = create_vector(4, {-1, -2, -3, -4});
    print(M);
    // y = 2.0 * M * x + 3.0 * y
    matrix_vector(2.0, M, false, x, 3.0, y);
    print(y);
    ASSERT_EQ(y, create_vector(4, {67, 214, 361, 508}));

    x = create_vector(5, {5, 4, 3, 2, 1});
    y = create_vector(4, {-1, -2, -3, -4});
    matrix_vector(2.0, M, true, y, 3.0, x);
    print(x);
    ASSERT_EQ(x, create_vector(5, {-205, -228, -251, -274, -297}));
}
