#include <iostream>
#include <vector>
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "matvec_oper.hpp"

TEST(test_matvec_oper, submatrix)
{
    matrix A = create_matrix(4, 3,
                             {
                                 1, 2, 3,  //
                                 2, 2, 1,  //
                                 7, 2, 9,  //
                                 7, 2, 9,  //
                             });
    auto A_c_m_v = matrix_col_const_view(A, 1);
    auto A_c_c_v = matrix_col_mutable_view(A, 1);
    auto A_r_m_v = matrix_row_mutable_view(A, 1);
    auto A_r_c_v = matrix_row_const_view(A, 1);
    print(A_c_m_v);
    print(A_c_c_v);
    print(A_r_m_v);
    print(A_r_c_v);
    ASSERT_EQ(A_c_m_v, create_vector(A.m(), {2, 2, 2, 2}));
    ASSERT_EQ(A_c_c_v, create_vector(A.m(), {2, 2, 2, 2}));
    ASSERT_EQ(A_r_m_v, create_vector(A.n(), {2, 2, 1}));
    ASSERT_EQ(A_r_c_v, create_vector(A.n(), {2, 2, 1}));
}
