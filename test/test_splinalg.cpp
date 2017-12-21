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
                                 CPR_COL);
        auto pi = vector(4, 1.0);
        auto b = vector(4, 0.0);
        real_t prec = linsv_sor_method(Q, true, pi, b, 1.0, 1.2, tol, max_iter,
                                       check_interval);
        print(prec);
        print(Q);
        print(pi);
        print(dot(Q, true, pi));
        ASSERT_TRUE(near_eq(dot(Q, true, pi), b, 1e-6));
    }
    {
        auto QTT = create_spmatrix(4, 4,
                                   {
                                       -2, 2, 0, 0,  //
                                       1, -2, 1, 0,  //
                                       0, 1, -2, 1,  //
                                       0, 0, 1, -2   //
                                   },
                                   CPR_COL);
        auto pi = create_vector(4.0, {1.0, 0, 0, 0});
        blas_scale(-1.0, pi);
        auto tau = vector(4, 0.0);
        linsv_sor_method(QTT, true, tau, pi, 1.2, tol, max_iter,
                         check_interval);
        print(tau);
        print(QTT);
        print(pi);
        ASSERT_TRUE(near_eq(dot(QTT, true, tau), pi, 1e-6));
    }
}
