#include "common.hpp"
#include "debug.hpp"
#include "gtest/gtest.h"
#include "linalg.hpp"
#include "matvec.hpp"

TEST(test_linalg, decomp_qr_gms1)
{
    auto A = create_matrix(3, 3,
                           {

                               12, -51, 4,   //
                               6, 167, -68,  //
                               -4, 24, -41,  //
                           });
    auto Q = A;
    auto R = matrix(A.n(), A.n());
    decomp_qr_mgs(Q, R);
    print(A);
    print(Q);
    print(R);
    ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
}
TEST(test_linalg, decomp_qr_gms2)
{
    auto A = create_matrix(3, 3,
                           {

                               12, -51, -12,  //
                               6, 167, -6,    //
                               -4, 24, 4,     //
                           });
    auto Q = A;
    auto R = matrix(A.n(), A.n());
    decomp_qr_mgs(Q, R);
    print(A);
    print(Q);
    print(R);
    ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
}
TEST(test_linalg, decomp_qr_gms3)
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
    decomp_qr_mgs(Q, R);
    print(A);
    print(Q);
    print(R);
    ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
}
TEST(test_linalg, householder1)
{
    auto w = create_vector(5, {1, 2, 3, 4, 5});
    auto v = w;
    real_t tau = find_householder_vector(v);
    apply_householder_reflector(w, tau, v);
    print(w);
    ASSERT_TRUE(
        near_eq(w, create_vector(5, {-blas_norm2(w), 0, 0, 0, 0}), 1e-6));
}
TEST(test_linalg, householder2)
{
    auto w = create_vector(5, {-1, 2, 3, 4, 5});
    auto v = w;
    real_t tau = find_householder_vector(v);
    apply_householder_reflector(w, tau, v);
    print(w);
    ASSERT_TRUE(
        near_eq(w, create_vector(5, {blas_norm2(w), 0, 0, 0, 0}), 1e-6));
}
TEST(test_linalg, householder3)
{
    auto A = create_matrix(3, 3,
                           {
                               12, -51, -12,  //
                               6, 167, -6,    //
                               -4, 24, 4,     //
                           });
    vector v(3);
    blas_copy(A.col(0), v);
    real_t tau = find_householder_vector(v);
    apply_householder_reflector(A, tau, v);
    print(A);
    ASSERT_TRUE(near_eq(A.col(0),
                        create_vector(3, {-blas_norm2(A.col(0)), 0, 0}), 1e-6));
}

TEST(test_linalg, decomp_qr_hr1)
{
    auto A = create_matrix(3, 3,
                           {

                               12, -51, 4,   //
                               6, 167, -68,  //
                               -4, 24, -41,  //
                           });
    auto QR = A;
    vector tau(A.n());
    decomp_qr_hr(QR, tau);
    matrix Q(A.m(), A.n());
    matrix R(A.n(), A.n());
    std::cout << "A ";
    print(A);
    std::cout << "QR ";
    print(QR);
    unpack_qr(QR, tau, Q, R);
    std::cout << "Q ";
    print(Q);
    std::cout << "R ";
    print(R);
    std::cout << "Q.R ";
    print(dot(Q, R));
    ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
}
TEST(test_linalg, decomp_qr_hr2)
{
    auto A = create_matrix(3, 3,
                           {

                               12, -51, -12,  //
                               6, 167, -6,    //
                               -4, 24, 4,     //
                           });
    auto QR = A;
    vector tau(A.n());
    decomp_qr_hr(QR, tau);
    matrix Q(A.m(), A.n());
    matrix R(A.n(), A.n());
    std::cout << "A ";
    print(A);
    std::cout << "QR ";
    print(QR);
    unpack_qr(QR, tau, Q, R);
    std::cout << "Q ";
    print(Q);
    std::cout << "R ";
    print(R);
    std::cout << "Q.R ";
    print(dot(Q, R));
    ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
}
TEST(test_linalg, decomp_qr_hr3)
{
    auto A = create_matrix(4, 3,
                           {

                               12, -51, 4,   //
                               6, 167, -68,  //
                               -4, 24, -41,  //
                               1, 2, 3       //
                           });
    auto QR = A;
    vector tau(A.n());
    decomp_qr_hr(QR, tau);
    matrix Q(A.m(), A.n());
    matrix R(A.n(), A.n());
    std::cout << "A ";
    print(A);
    std::cout << "QR ";
    print(QR);
    unpack_qr(QR, tau, Q, R);
    std::cout << "Q ";
    print(Q);
    std::cout << "R ";
    print(R);
    std::cout << "Q.R ";
    print(dot(Q, R));
    ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
}

TEST(test_linalg, qr_compact_oper)
{
    auto A = create_matrix(4, 3,
                           {

                               12, -51, 4,   //
                               6, 167, -68,  //
                               -4, 24, -41,  //
                               -2, 5, 8      //
                           });
    const auto v = create_vector(4, {1, 2, 3, 4});
    auto QR = A;
    vector tau(A.n());
    decomp_qr_hr(QR, tau);
    matrix Q(A.m(), A.n());
    matrix R(A.n(), A.n());
    unpack_qr(QR, tau, Q, R);
    auto Qtv = v;
    qt_dot_vector(QR, tau, Qtv);
    std::cout << "compact ";
    print(Qtv.sub(0, 3));
    std::cout << "direct ";
    print(dot(Q.transpose(), v));
    ASSERT_TRUE(near_eq(Qtv.sub(0, 3), dot(Q.transpose(), v), 1e-6));

    const auto w = create_vector(4, {1, 2, 3, 0});
    auto Qw = w;
    q_dot_vector(QR, tau, Qw);
    std::cout << "compact ";
    print(Qw);
    std::cout << "direct ";
    print(dot(Q, w.sub(0, 3)));
    ASSERT_TRUE(near_eq(Qw, dot(Q, w.sub(0, 3)), 1e-6));
}

TEST(test_linalg, solve_upper_tri)
{
    auto U = create_matrix(4, 4,
                           {
                               1, 2, 0, 4,  //
                               0, 6, 8, 0,  //
                               0, 0, 2, 3,  //
                               0, 0, 0, 5   //
                           });
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    solve_upper_tri(U, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(U, x), b, 1e-6));
}
TEST(test_linalg, solve_lower_tri)
{
    auto U = create_matrix(4, 4,
                           {
                               1, 2, 0, 4,  //
                               0, 6, 8, 0,  //
                               0, 0, 2, 3,  //
                               0, 0, 0, 5   //
                           });
    auto L = U.transpose();
    const vector b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    solve_lower_tri(L, x);
    print(x);
    ASSERT_TRUE(near_eq(dot(L, x), b, 1e-6));
}

TEST(test_linalg, least_square)
{
    const auto A = create_matrix(4, 3,
                                 {

                                     12, -51, 4,   //
                                     6, 167, -68,  //
                                     -4, 24, -41,  //
                                     -2, 5, 8      //
                                 });
    const auto b = create_vector(4, {1, 2, 3, 4});
    vector x = b;
    matrix QR = A;
    least_square_qr(QR, x);
    print(x.sub(0, A.n()));
    ASSERT_TRUE(near_eq(x.sub(0, A.n()),
                        create_vector(3, {-0.0243225, -0.00982752, -0.0549788}),
                        1e-6));
}

TEST(test_linalg, decomp_lup)
{
    const auto A = create_matrix(3, 3,
                                 {

                                     12, -51, 4,   //
                                     6, 167, -68,  //
                                     -4, 24, -41,  //
                                 });
    matrix P = identity_matrix(3);
    matrix U = A;
    decomp_lup(U, P);
    matrix L(A.m(), A.n());
    unpack_lu(U, L);
    ASSERT_TRUE(near_eq(P * A, L * U, 1e-6));
}

TEST(test_linalg, decomp_lu)
{
    const auto A = create_matrix(3, 3,
                                 {

                                     12, -51, 4,   //
                                     6, 167, -68,  //
                                     -4, 24, -41,  //
                                 });
    matrix P = identity_matrix(3);
    matrix U = A;
    decomp_lu(U);
    matrix L(A.m(), A.n());
    unpack_lu(U, L);
    ASSERT_TRUE(near_eq(A, L * U, 1e-6));
}

TEST(test_linalg, solve_lu)
{
    const auto A_ = create_matrix(3, 3,
                                  {

                                      12, -51, 4,   //
                                      6, 167, -68,  //
                                      -4, 24, -41,  //
                                  });
    const auto b_ = create_vector(3, {1, 2, 3});
    matrix A = A_;
    vector x = b_;
    solve_lu(A, x);
    ASSERT_TRUE(near_eq(A_ * x, b_, 1e-6));
}

TEST(test_linalg, solve_qr)
{
    const auto A_ = create_matrix(3, 3,
                                  {

                                      12, -51, 4,   //
                                      6, 167, -68,  //
                                      -4, 24, -41,  //
                                  });
    const auto b_ = create_vector(3, {1, 2, 3});
    matrix A = A_;
    vector x = b_;
    solve_qr(A, x);
    ASSERT_TRUE(near_eq(A_ * x, b_, 1e-6));
}
