#include <iostream>
#include <vector>
#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "matvec.hpp"
TEST(test_matrix, submatrix)
{
    {
        matrix A = create_matrix(4, 3,
                                 {
                                     1, 2, 3,  //
                                     2, 2, 1,  //
                                     7, 2, 9,  //
                                     7, 2, 9,  //
                                 });
        auto A_c = A.col(1).sub(2);
        auto A_r = A.row(1);
        print(A_c);
        print(A_r);
        ASSERT_EQ(A_c, create_vector(A.m() - 2, {2, 2}));
        ASSERT_EQ(A_r, create_vector(A.n(), {2, 2, 1}));
        auto subA = A.sub(1, 1);
        ASSERT_EQ(subA, create_matrix(3, 2,
                                      {
                                          2, 1,  //
                                          2, 9,  //
                                          2, 9,  //
                                      }));
    }
    {
        matrix A_ = create_matrix(4, 3,
                                  {
                                      1, 2, 3,  //
                                      2, 2, 1,  //
                                      7, 2, 9,  //
                                      7, 2, 9,  //
                                  });
        auto A = matrix_mutable_view(A_);
        auto A_c = A.col(1).sub(2);
        auto A_r = A.row(1);
        print(A_c);
        print(A_r);
        ASSERT_EQ(A_c, create_vector(A.m() - 2, {2, 2}));
        ASSERT_EQ(A_r, create_vector(A.n(), {2, 2, 1}));
        auto subA = A.sub(1, 1);
        ASSERT_EQ(subA, create_matrix(3, 2,
                                      {
                                          2, 1,  //
                                          2, 9,  //
                                          2, 9,  //
                                      }));
    }
    {
        matrix A_ = create_matrix(4, 3,
                                  {
                                      1, 2, 3,  //
                                      2, 2, 1,  //
                                      7, 2, 9,  //
                                      7, 2, 9,  //
                                  });
        copy(create_vector(2, {8, 8}), A_.col(2).sub(2));
        ASSERT_EQ(A_.col(2), create_vector(4, {3, 1, 8, 8}));
    }
}
