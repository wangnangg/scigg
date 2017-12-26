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
    const matrix M = create_matrix(4, 5,
                                   {
                                       1,  2,  3,  4,  5,   //
                                       6,  7,  8,  9,  10,  //
                                       11, 12, 13, 14, 15,  //
                                       16, 17, 18, 19, 20   //
                                   });
    const matrix MT = create_matrix(5, 4,
                                    {
                                        1, 6,  11, 16,  //
                                        2, 7,  12, 17,  //
                                        3, 8,  13, 18,  //
                                        4, 9,  14, 19,  //
                                        5, 10, 15, 20   //
                                    });
    ASSERT_EQ(M.transpose(), MT);

    {
        vector x = create_vector(5, {5, 4, 3, 2, 1});
        vector y = create_vector(4, {-1, -2, -3, -4});
        // y = 2.0 * M * x + 3.0 * y
        blas_matrix_vector(2.0, M, x, 3.0, y);
        ASSERT_EQ(y, create_vector(4, {67, 214, 361, 508}));
    }
    {
        vector x = create_vector(5, {5, 4, 3, 2, 1});
        vector y = create_vector(4, {-1, -2, -3, -4});
        // y = 2.0 * M * x + 3.0 * y
        blas_matrix_vector(2.0, MT.transpose(), x, 3.0, y);
        ASSERT_EQ(y, create_vector(4, {67, 214, 361, 508}));
    }
    {
        vector x = create_vector(5, {5, 4, 3, 2, 1});
        vector y = create_vector(4, {-1, -2, -3, -4});
        blas_matrix_vector(2.0, M.transpose(), y, 3.0, x);
        ASSERT_EQ(x, create_vector(5, {-205, -228, -251, -274, -297}));
    }
    {
        vector x = create_vector(5, {5, 4, 3, 2, 1});
        vector y = create_vector(4, {-1, -2, -3, -4});
        blas_matrix_vector(2.0, MT, y, 3.0, x);
        ASSERT_EQ(x, create_vector(5, {-205, -228, -251, -274, -297}));
    }
}

TEST(test_blas, level3)
{
    const matrix A_ = create_matrix(4, 5,
                                    {
                                        1,  2,  3,  4,  5,   //
                                        6,  7,  8,  9,  10,  //
                                        11, 12, 13, 14, 15,  //
                                        16, 17, 18, 19, 20   //
                                    });

    const matrix AT_ = create_matrix(5, 4,
                                     {
                                         1, 6,  11, 16,  //
                                         2, 7,  12, 17,  //
                                         3, 8,  13, 18,  //
                                         4, 9,  14, 19,  //
                                         5, 10, 15, 20   //
                                     });
    ASSERT_EQ(A_, AT_.transpose());
    ASSERT_EQ(A_.transpose(), AT_);
    const matrix B_ = create_matrix(5, 3,
                                    {
                                        1, -2, 3,   //
                                        6, -6, 8,   //
                                        11, 5, 13,  //
                                        16, 4, 18,  //
                                        1, -4, -3,  //
                                    });
    const matrix BT_ = create_matrix(3, 5,
                                     {
                                         1, 6, 11, 16, 1,   //
                                         -2, -6, 5, 4, -4,  //
                                         3, 8, 13, 18, -3   //
                                     });
    ASSERT_EQ(B_, BT_.transpose());
    ASSERT_EQ(B_.transpose(), BT_);
    const matrix C_ = create_matrix(4, 3,
                                    {
                                        1, 1, 1,   //
                                        2, 2, 2,   //
                                        3, 5, 13,  //
                                        4, 4, 18,  //
                                    });
    const matrix CT_ = create_matrix(3, 4,
                                     {
                                         1, 2, 3, 4,    //
                                         1, 2, 5, 4,    //
                                         1, 2, 13, 18,  //
                                     });
    ASSERT_EQ(C_, CT_.transpose());
    ASSERT_EQ(C_.transpose(), CT_);

    {
        // C = alpha * A * B + beta * C
        matrix C = C_;
        blas_matrix_matrix(2.0, A_, B_, 3.0, C);
        ASSERT_EQ(C, create_matrix(4, 3,
                                   {
                                       233, -3, 233,     //
                                       586, -30, 626,    //
                                       939, -51, 1049,   //
                                       1292, -84., 1454  //
                                   }));
    }
    {
        // C = alpha * A * B + beta * C
        matrix C = C_;
        blas_matrix_matrix(2.0, AT_.transpose(), B_, 3.0, C);
        ASSERT_EQ(C, create_matrix(4, 3,
                                   {
                                       233, -3, 233,     //
                                       586, -30, 626,    //
                                       939, -51, 1049,   //
                                       1292, -84., 1454  //
                                   }));
    }
    {
        // C = alpha * A * B + beta * C
        matrix C = C_;
        blas_matrix_matrix(2.0, A_, BT_.transpose(), 3.0, C);
        ASSERT_EQ(C, create_matrix(4, 3,
                                   {
                                       233, -3, 233,     //
                                       586, -30, 626,    //
                                       939, -51, 1049,   //
                                       1292, -84., 1454  //
                                   }));
    }
    {
        // C = alpha * A * B + beta * C
        matrix C = C_;
        blas_matrix_matrix(2.0, AT_.transpose(), BT_.transpose(), 3.0, C);
        ASSERT_EQ(C, create_matrix(4, 3,
                                   {
                                       233, -3, 233,     //
                                       586, -30, 626,    //
                                       939, -51, 1049,   //
                                       1292, -84., 1454  //
                                   }));
    }
    {
        // C = alpha * A * B + beta * C
        matrix C = CT_;
        blas_matrix_matrix(2.0, A_, B_, 3.0, C.transpose());
        ASSERT_EQ(C.transpose(), create_matrix(4, 3,
                                               {
                                                   233, -3, 233,     //
                                                   586, -30, 626,    //
                                                   939, -51, 1049,   //
                                                   1292, -84., 1454  //
                                               }));
    }
    {
        // C = alpha * A * B + beta * C
        matrix C = CT_;
        blas_matrix_matrix(2.0, AT_.transpose(), B_, 3.0, C.transpose());
        ASSERT_EQ(C.transpose(), create_matrix(4, 3,
                                               {
                                                   233, -3, 233,     //
                                                   586, -30, 626,    //
                                                   939, -51, 1049,   //
                                                   1292, -84., 1454  //
                                               }));
    }
    {
        // C = alpha * A * B + beta * C
        matrix C = CT_;
        blas_matrix_matrix(2.0, A_, BT_.transpose(), 3.0, C.transpose());
        ASSERT_EQ(C.transpose(), create_matrix(4, 3,
                                               {
                                                   233, -3, 233,     //
                                                   586, -30, 626,    //
                                                   939, -51, 1049,   //
                                                   1292, -84., 1454  //
                                               }));
    }
    {
        // C = alpha * A * B + beta * C
        matrix C = CT_;
        blas_matrix_matrix(2.0, AT_.transpose(), BT_.transpose(), 3.0,
                           C.transpose());
        ASSERT_EQ(C.transpose(), create_matrix(4, 3,
                                               {
                                                   233, -3, 233,     //
                                                   586, -30, 626,    //
                                                   939, -51, 1049,   //
                                                   1292, -84., 1454  //
                                               }));
    }
}
