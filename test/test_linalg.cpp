#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "linalg.hpp"
#include "matvec_oper.hpp"

TEST(test_linalg, qr_decomp_gms)
{
    {
        auto A = create_matrix(3, 3,
                               {

                                   12, -51, 4,   //
                                   6, 167, -68,  //
                                   -4, 24, -41,  //
                               });
        auto Q = A;
        auto R = matrix(A.n(), A.n());
        qr_decomp_mgs(Q, R);
        print(A);
        print(Q);
        print(R);
        ASSERT_TRUE(near_eq(A, dot(Q, false, R, false), 1e-8));
    }
    {
        auto A = create_matrix(3, 3,
                               {

                                   12, -51, -12,  //
                                   6, 167, -6,    //
                                   -4, 24, 4,     //
                               });
        auto Q = A;
        auto R = matrix(A.n(), A.n());
        qr_decomp_mgs(Q, R);
        print(A);
        print(Q);
        print(R);
        ASSERT_TRUE(near_eq(A, dot(Q, false, R, false), 1e-8));
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
        auto R = matrix(A.n(), A.n());
        qr_decomp_mgs(Q, R);
        print(A);
        print(Q);
        print(R);
        ASSERT_TRUE(near_eq(A, dot(Q, false, R, false), 1e-8));
    }
}

TEST(test_linalg, qr_decomp_hr)
{
    {
        auto A = create_matrix(3, 3,
                               {

                                   12, -51, 4,   //
                                   6, 167, -68,  //
                                   -4, 24, -41,  //
                               });
        auto R = A;
        auto V = matrix(A.m(), A.n());
        qr_decomp_hr(R, V);
        print(A);
        print(R);
        print(V);
        matrix Q(A.m(), A.m());
        recover_q_from_v(V, Q);
        print(Q);
        ASSERT_TRUE(near_eq(A, dot(Q, false, R, false), 1e-8));
    }
    {
        auto A = create_matrix(3, 3,
                               {

                                   12, -51, 4,   //
                                   6, 167, -68,  //
                                   -4, 24, -41,  //
                               });
        auto R = A;
        auto V = matrix(A.m(), A.n());
        qr_decomp_hr(R, V);
        print(A);
        print(R);
        print(V);
        matrix Q(A.m(), A.m());
        recover_q_from_v(V, Q);
        print(Q);
        ASSERT_TRUE(near_eq(A, dot(Q, false, R, false), 1e-8));
    }
    {
        auto A = create_matrix(4, 3,
                               {

                                   12, -51, 4,   //
                                   6, 167, -68,  //
                                   -4, 24, -41,  //
                                   1, 2, 3       //
                               });
        auto R = A;
        auto V = matrix(A.m(), A.n());
        qr_decomp_hr(R, V);
        print(A);
        print(R);
        print(V);
        matrix Q(A.m(), A.m());
        recover_q_from_v(V, Q);
        print(Q);
        ASSERT_TRUE(near_eq(A, dot(Q, false, R, false), 1e-8));
    }
}
