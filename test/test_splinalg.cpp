#include "blas.hpp"
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "linalg.hpp"
#include "matvec_oper.hpp"
#include "splinalg.hpp"
#include "spmatrix_oper.hpp"

static real_t tol = 1e-10;
static uint_t max_iter = 1000;
static uint_t check_interval = 10;

TEST(test_splinalg, sor_method_sum)
{
    auto Q = create_spmatrix(4, 4,
                             {
                                 -2, 2, 0, 0,  //
                                 1, -3, 2, 0,  //
                                 0, 1, -3, 2,  //
                                 1, 1, 1, 1    //
                             },
                             false);
    auto pi = vector(4, 0.0);
    auto b = vector(4, 0.0);
    b[3] = 1.0;
    real_t prec = spsolve_sor_method(Q.transpose(), pi, b, 0.8, tol, max_iter,
                                     check_interval);
    print(Q);
    print(pi);
    print(dot(Q.transpose(), pi));
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(Q.transpose(), pi), b, 1e-6));
}

TEST(test_splinalg, sor_method1)
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
    auto x = vector(pi.dim(), 0.0);
    real_t w = 1.2;
    real_t prec = spsolve_sor_method(QTT.transpose(), x, pi, w, tol, max_iter,
                                     check_interval);
    ASSERT_LT(prec, tol);
    print(prec);
    print(x);
    ASSERT_TRUE(near_eq(dot(QTT.transpose(), x), pi, 1e-2));
}
TEST(test_splinalg, sor_method2)
{
    auto QTT = create_spmatrix(4, 4,
                               {
                                   -4, 2, 0, 1,  //
                                   1, -6, 1, 3,  //
                                   2, 1, -5, 1,  //
                                   3, 2, 1, -6   //
                               },
                               false);
    auto pi = create_vector(4.0, {1.0, 0, 0, 0});
    blas_scale(-1.0, pi);
    auto x = vector(pi.dim(), 0.0);
    real_t w = 1.2;
    real_t prec = spsolve_sor_method(QTT.transpose(), x, pi, w, tol, max_iter,
                                     check_interval);
    ASSERT_LT(prec, tol);
    print(prec);
    print(x);
    ASSERT_TRUE(near_eq(dot(QTT.transpose(), x), pi, 1e-2));
}

TEST(test_splinalg, power_method1)
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
    real_t prec =
        eigen_power_method(P.transpose(), pi, tol, max_iter, check_interval);
    print(P);
    print(pi);
    print(dot(P.transpose(), pi));
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(P.transpose(), pi), pi, 1e-6));
}
TEST(test_splinalg, power_method2)
{
    auto P = create_spmatrix(4, 4,
                             {
                                 0.8, 0.2, 0, 0,    //
                                 0.3, 0.5, 0.2, 0,  //
                                 0, 0.3, 0.5, 0.2,  //
                                 0, 0, 0.3, 0.7     //
                             },
                             true);
    auto pi = vector(4, 0.0);
    pi[0] = 1.0;
    real_t prec =
        eigen_power_method(P.transpose(), pi, tol, max_iter, check_interval);
    print(P);
    print(pi);
    print(dot(P.transpose(), pi));
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(P.transpose(), pi), pi, 1e-6));
}

TEST(test_splinalg, gmres_method)
{
    auto QTT = create_spmatrix(4, 4,
                               {
                                   -4, 2, 0, 1,  //
                                   1, -6, 1, 3,  //
                                   2, 1, -5, 1,  //
                                   3, 2, 1, -6   //
                               },
                               false);
    auto pi = create_vector(4.0, {1.0, 0, 0, 0});
    scale(pi, -1.0);
    auto tau = vector(4, 0.0);
    spsolve_gmres_gms(QTT.transpose(), tau, pi, 4, 1e-6, 4);

    std::cout << "tau ";
    print(tau);
    print(QTT);
    print(pi);
    ASSERT_TRUE(near_eq(dot(QTT.transpose(), tau), pi, 1e-6));
}

TEST(test_splinalg, restart_gmres_method)
{
    auto QTT = create_spmatrix(4, 4,
                               {
                                   -4, 2, 0, 1,  //
                                   1, -6, 1, 3,  //
                                   2, 1, -5, 1,  //
                                   3, 2, 1, -6   //
                               },
                               false);
    auto pi = create_vector(4.0, {1.0, 0, 0, 0});
    scale(pi, -1.0);
    auto tau = vector(4, 0.0);
    real_t prec =
        spsolve_restart_gmres_gms(QTT.transpose(), tau, pi, 2, tol, 100, 2);
    std::cout << "tau ";
    print(tau);
    print(QTT);
    print(pi);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(QTT.transpose(), tau), pi, tol));
}

TEST(test_splinalg, spdecomp_ilu1)
{
    const auto A = create_spmatrix(4, 4,
                                   {
                                       -4, 2, 0, 1,  //
                                       1, -6, 1, 3,  //
                                       2, 1, -5, 1,  //
                                       3, 2, 1, -6   //
                                   },
                                   true);
    auto iLU = A;
    spdecomp_ilu(iLU);
    auto LU = spmatrix2dense(A);
    decomp_lu(LU);
    ASSERT_TRUE(near_eq(spmatrix2dense(iLU), LU, 1e-6));
}

TEST(test_splinalg, spdecomp_ilu2)
{
    const auto A = create_spmatrix(4, 4,
                                   {
                                       1, 2, 0, 1,   //
                                       0, -6, 0, 3,  //
                                       2, 0, -5, 0,  //
                                       1, 0, 1, -6   //
                                   },
                                   true);
    auto iLU = A;
    spdecomp_ilu(iLU);
    auto LU = spmatrix2dense(A);
    decomp_lu(LU);
    for (size_t i = 0; i < iLU.m(); i++)
    {
        auto view = iLU[i];
        for (size_t j = 0; j < view.nnz; j++)
        {
            near_eq(LU(i, view.idx[j]), view.val[j], 1e-6);
        }
    }
}
TEST(test_splinalg, spsolve_lower_tri1)
{
    auto L = create_spmatrix(4, 4,
                             {
                                 1, 0, 0, 0,  //
                                 2, 6, 0, 0,  //
                                 0, 8, 2, 0,  //
                                 4, 0, 3, 5   //
                             },
                             true);
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_lower_tri(L, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(L, x), b, 1e-6));
}

TEST(test_splinalg, spsolve_lower_tri2)
{
    auto L = create_spmatrix(4, 4,
                             {
                                 1, 0, 0, 0,  //
                                 2, 6, 0, 0,  //
                                 0, 8, 2, 0,  //
                                 4, 0, 3, 5   //
                             },
                             false);
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_lower_tri(L, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(L, x), b, 1e-6));
}

TEST(test_splinalg, spsolve_lower_tri3)
{
    auto U = create_spmatrix(4, 4,
                             {
                                 1, 2, 0, 4,  //
                                 0, 6, 8, 0,  //
                                 0, 0, 2, 3,  //
                                 0, 0, 0, 5   //
                             },
                             true);
    auto L = U.transpose();
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_lower_tri(L, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(L, x), b, 1e-6));
}

TEST(test_splinalg, spsolve_lower_tri4)
{
    auto U = create_spmatrix(4, 4,
                             {
                                 1, 2, 0, 4,  //
                                 0, 6, 8, 0,  //
                                 0, 0, 2, 3,  //
                                 0, 0, 0, 5   //
                             },
                             false);
    auto L = U.transpose();
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_lower_tri(L, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(L, x), b, 1e-6));
}

TEST(test_splinalg, spsolve_upper_tri1)
{
    auto L = create_spmatrix(4, 4,
                             {
                                 1, 0, 0, 0,  //
                                 2, 6, 0, 0,  //
                                 0, 8, 2, 0,  //
                                 4, 0, 3, 5   //
                             },
                             true);
    auto U = L.transpose();
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_upper_tri(U, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(U, x), b, 1e-6));
}

TEST(test_splinalg, spsolve_upper_tri2)
{
    auto L = create_spmatrix(4, 4,
                             {
                                 1, 0, 0, 0,  //
                                 2, 6, 0, 0,  //
                                 0, 8, 2, 0,  //
                                 4, 0, 3, 5   //
                             },
                             false);
    auto U = L.transpose();
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_upper_tri(U, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(U, x), b, 1e-6));
}

TEST(test_splinalg, spsolve_upper_tri3)
{
    auto U = create_spmatrix(4, 4,
                             {
                                 1, 2, 0, 4,  //
                                 0, 6, 8, 0,  //
                                 0, 0, 2, 3,  //
                                 0, 0, 0, 5   //
                             },
                             true);
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_upper_tri(U, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(U, x), b, 1e-6));
}

TEST(test_splinalg, spsolve_upper_tri4)
{
    auto U = create_spmatrix(4, 4,
                             {
                                 1, 2, 0, 4,  //
                                 0, 6, 8, 0,  //
                                 0, 0, 2, 3,  //
                                 0, 0, 0, 5   //
                             },
                             false);
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    spsolve_upper_tri(U, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(U, x), b, 1e-6));
}
