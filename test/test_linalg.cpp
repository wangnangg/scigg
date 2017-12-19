#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "linalg.hpp"
#include "matvec_oper.hpp"

TEST(test_linalg, QR_decomp_GMS)
{
    {
        auto A = create_matrix(3, 3,
                               {

                                   12, -51, 4,   //
                                   6, 167, -68,  //
                                   -4, 24, -41,  //
                               });
        auto Q = A;
        auto R = matrix(3, 3);
        QR_decomp_MGS(Q, R);
        print(A);
        print(Q);
        print(R);
        ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
    }
    {
        auto A = create_matrix(3, 3,
                               {

                                   12, -51, -12,  //
                                   6, 167, -6,    //
                                   -4, 24, 4,     //
                               });
        auto Q = A;
        auto R = matrix(3, 3);
        QR_decomp_MGS(Q, R);
        print(A);
        print(Q);
        print(R);
        ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
    }
    {
        auto A = create_matrix(4, 3,
                               {

                                   12, -51, 4,   //
                                   6, 167, -68,  //
                                   -4, 24, -41,  //
                                   1, 2, 3       //
                               });
        auto Q = A;
        auto R = matrix(3, 3);
        QR_decomp_MGS(Q, R);
        print(A);
        print(Q);
        print(R);
        ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
    }
}
