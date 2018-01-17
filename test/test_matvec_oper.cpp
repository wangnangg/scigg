#include <iostream>
#include <vector>
#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "matvec.hpp"

TEST(test_matvec_oper, vector_oper)
{
    auto x = create_vector(4, {1, 2, 3, 4});
    auto y = create_vector(4, {2, 3, 4, 5});
    ASSERT_EQ(x + y, create_vector(4, {3, 5, 7, 9}));
    ASSERT_EQ(x - y, create_vector(4, {-1, -1, -1, -1}));
    ASSERT_EQ(x * y, 2 + 6 + 12 + 20);
    x += y;
    ASSERT_EQ(x, create_vector(4, {3, 5, 7, 9}));
    x -= y;
    ASSERT_EQ(x, create_vector(4, {1, 2, 3, 4}));
}
TEST(test_matvec_oper, matrix_oper)
{
    const matrix A = create_matrix(4, 5,
                                   {
                                       1,  2,  3,  4,  5,   //
                                       6,  7,  8,  9,  10,  //
                                       11, 12, 13, 14, 15,  //
                                       16, 17, 18, 19, 20   //
                                   });
    const matrix B = create_matrix(5, 3,
                                   {
                                       1, -2, 3,   //
                                       6, -6, 8,   //
                                       11, 5, 13,  //
                                       16, 4, 18,  //
                                       1, -4, -3,  //
                                   });
    ASSERT_EQ(A * B, create_matrix(4, 3,
                                   {115, -3, 115, 290, -18, 310, 465, -33, 505,
                                    640, -48, 700}));
}
TEST(teset_matvec_oper, matrix_vector_oper)
{
    const matrix A = create_matrix(4, 5,
                                   {
                                       1,  2,  3,  4,  5,   //
                                       6,  7,  8,  9,  10,  //
                                       11, 12, 13, 14, 15,  //
                                       16, 17, 18, 19, 20   //
                                   });
    auto x = create_vector(5, {1, 2, 3, 4, 5});
    ASSERT_EQ(A * x, create_vector(4, {55, 130, 205, 280}));
}
