#include <iostream>
#include <vector>
#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "matvec.hpp"
#include "spmatvec.hpp"

TEST(test_spblas, level2_crs)
{
    spmatrix M = create_spmatrix(4, 5,
                                 {
                                     1,  2,  3,  4,  5,   //
                                     6,  7,  8,  9,  10,  //
                                     11, 12, 13, 14, 15,  //
                                     16, 17, 18, 19, 20   //
                                 },
                                 true);

    vector x = create_vector(5, {5, 4, 3, 2, 1});
    vector y = create_vector(4, {-1, -2, -3, -4});
    print(M);
    // y = 2.0 * M * x + 3.0 * y
    spblas_matrix_vector(2.0, M, x, 3.0, y);
    print(y);
    ASSERT_EQ(y, create_vector(4, {67, 214, 361, 508}));

    x = create_vector(5, {5, 4, 3, 2, 1});
    y = create_vector(4, {-1, -2, -3, -4});
    spblas_matrix_vector(2.0, M.transpose(), y, 3.0, x);
    print(x);
    ASSERT_EQ(x, create_vector(5, {-205, -228, -251, -274, -297}));
}
TEST(test_spblas, level2_ccs)
{
    spmatrix M = create_spmatrix(4, 5,
                                 {
                                     1,  2,  3,  4,  5,   //
                                     6,  7,  8,  9,  10,  //
                                     11, 12, 13, 14, 15,  //
                                     16, 17, 18, 19, 20   //
                                 },
                                 false);

    vector x = create_vector(5, {5, 4, 3, 2, 1});
    vector y = create_vector(4, {-1, -2, -3, -4});
    print(M);
    // y = 2.0 * M * x + 3.0 * y
    spblas_matrix_vector(2.0, M, x, 3.0, y);
    print(y);
    ASSERT_EQ(y, create_vector(4, {67, 214, 361, 508}));

    x = create_vector(5, {5, 4, 3, 2, 1});
    y = create_vector(4, {-1, -2, -3, -4});
    spblas_matrix_vector(2.0, M.transpose(), y, 3.0, x);
    print(x);
    ASSERT_EQ(x, create_vector(5, {-205, -228, -251, -274, -297}));
}
