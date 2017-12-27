#include "blas.hpp"
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "matvec_oper.hpp"
#include "splinalg.hpp"
#include "spmatrix_oper.hpp"

static real_t tol = 1e-10;
static int_t max_iter = 1000;
static int_t check_interval = 10;
TEST(test_splinalg, sor_method)
{
    {
        auto Q = create_spmatrix(4, 4,
                                 {
                                     -2, 2, 0, 0,  //
                                     1, -3, 2, 0,  //
                                     0, 1, -3, 2,  //
                                     0, 0, 1, -1   //
                                 },
                                 false);
        auto pi = vector(4, 1.0);
        auto b = vector(4, 0.0);
        real_t prec = linsv_sor_method(Q.transpose(), pi, b, 1.0, 1.2, tol,
                                       max_iter, check_interval);
        print(Q);
        print(pi);
        print(dot(Q.transpose(), pi));
        ASSERT_LT(prec, tol);
        ASSERT_TRUE(near_eq(dot(Q.transpose(), pi), b, 1e-6));
    }
    {
        auto QTT = create_spmatrix(4, 4,
                                   {
                                       -2, 2, 0, 0,  //
                                       1, -2, 1, 0,  //
                                       0, 1, -2, 1,  //
                                       0, 0, 1, -2   //
                                   },
                                   false);
        auto pi = create_vector(4.0, {1.0, 0, 0, 0});
        blas_scale(-1.0, pi);
        auto tau = vector(4, 0.0);
        real_t prec = linsv_sor_method(QTT.transpose(), tau, pi, 1.2, tol,
                                       max_iter, check_interval);
        print(tau);
        print(QTT);
        print(pi);
        ASSERT_LT(prec, tol);
        ASSERT_TRUE(near_eq(dot(QTT.transpose(), tau), pi, 1e-6));
    }
}

TEST(test_splinalg, power_method)
{
    {
        auto P = create_spmatrix(4, 4,
                                 {
                                     0.8, 0.2, 0, 0,    //
                                     0.3, 0.5, 0.2, 0,  //
                                     0, 0.3, 0.5, 0.2,  //
                                     0, 0, 0.3, 0.7     //
                                 },
                                 false);
        auto pi = vector(4, 0.0);
        pi[0] = 1.0;
        real_t prec = eigen_power_method(P.transpose(), pi, tol, max_iter,
                                         check_interval);
        print(P);
        print(pi);
        print(dot(P.transpose(), pi));
        ASSERT_LT(prec, tol);
        ASSERT_TRUE(near_eq(dot(P.transpose(), pi), pi, 1e-6));
    }
    {
        auto P = create_spmatrix(4, 4,
                                 {
                                     0.8, 0.2, 0, 0,    //
                                     0.3, 0.5, 0.2, 0,  //
                                     0, 0.3, 0.5, 0.2,  //
                                     0, 0, 0.3, 0.7     //
                                 },
                                 false);
        auto pi = vector(4, 0.0);
        pi[0] = 1.0;
        real_t prec = eigen_power_method(P.transpose(), pi, tol, max_iter,
                                         check_interval);
        print(P);
        print(pi);
        print(dot(P.transpose(), pi));
        ASSERT_LT(prec, tol);
        ASSERT_TRUE(near_eq(dot(P.transpose(), pi), pi, 1e-6));
    }
}
