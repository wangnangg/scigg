#include <iostream>
#include <vector>
#include "blas.hpp"
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "matvec_oper.hpp"

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

    ASSERT_EQ(blas_dot(v1, v2), 32);

    vector v3(2, 0);
    v3[0] = 3;
    v3[1] = 4;
    ASSERT_EQ(blas_norm2(v3), 5);

    ASSERT_EQ(blas_abs_sum(v1), 6);

    ASSERT_EQ(blas_abs_max_idx(v2), 2);

    blas_swap(v1, v2);
    ASSERT_EQ(v1[0], 4);
    ASSERT_EQ(v1[1], 5);
    ASSERT_EQ(v1[2], 6);
    blas_swap(v1, v2);

    vector v4(3, 0);
    blas_copy(v1, v4);
    ASSERT_TRUE(v4 == v1);

    blas_axpy(2.0, v1, v2);
    ASSERT_EQ(v2[0], 6);
    ASSERT_EQ(v2[1], 9);
    ASSERT_EQ(v2[2], 12);
    blas_axpy(-2.0, v1, v2);

    blas_scale(10, v1);
    ASSERT_EQ(v1[0], 10);
    ASSERT_EQ(v1[1], 20);
    ASSERT_EQ(v1[2], 30);
    blas_scale(0.1, v1);
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
    blas_matrix_vector(2.0, M, false, x, 3.0, y);
    print(y);
    ASSERT_EQ(y, create_vector(4, {67, 214, 361, 508}));

    x = create_vector(5, {5, 4, 3, 2, 1});
    y = create_vector(4, {-1, -2, -3, -4});
    blas_matrix_vector(2.0, M, true, y, 3.0, x);
    print(x);
    ASSERT_EQ(x, create_vector(5, {-205, -228, -251, -274, -297}));
}

TEST(test_blas, level3)
{
    matrix A = create_matrix(4, 5,
                             {
                                 1,  2,  3,  4,  5,   //
                                 6,  7,  8,  9,  10,  //
                                 11, 12, 13, 14, 15,  //
                                 16, 17, 18, 19, 20   //
                             });
    matrix B = create_matrix(5, 3,
                             {
                                 1, -2, 3,   //
                                 6, -6, 8,   //
                                 11, 5, 13,  //
                                 16, 4, 18,  //
                                 1, -4, -3,  //
                             });
    matrix C = create_matrix(4, 3,
                             {
                                 1, 1, 1,   //
                                 2, 2, 2,   //
                                 3, 5, 13,  //
                                 4, 4, 18,  //
                             });

    // C = alpha * A * B + beta * C
    blas_matrix_matrix(2.0, A, false, B, false, 3.0, C);
    ASSERT_EQ(C, create_matrix(4, 3,
                               {
                                   233, -3, 233,     //
                                   586, -30, 626,    //
                                   939, -51, 1049,   //
                                   1292, -84., 1454  //
                               }));
}