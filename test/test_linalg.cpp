#include "blas.hpp"
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
        auto R = matrix(A.n(), A.n());
        qr_decomp_mgs(Q, R);
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
        auto R = matrix(A.n(), A.n());
        qr_decomp_mgs(Q, R);
        print(A);
        print(Q);
        print(R);
        ASSERT_TRUE(near_eq(A, dot(Q, R), 1e-8));
    }
}
TEST(test_linalg, householder)
{
    {
        auto w = create_vector(5, {1, 2, 3, 4, 5});
        auto v = w;
        real_t tau = find_householder_vector(v);
        apply_householder_reflector(w, tau, v);
        print(w);
        ASSERT_TRUE(
            near_eq(w, create_vector(5, {-blas_norm2(w), 0, 0, 0, 0}), 1e-6));
    }
    {
        auto w = create_vector(5, {-1, 2, 3, 4, 5});
        auto v = w;
        real_t tau = find_householder_vector(v);
        apply_householder_reflector(w, tau, v);
        print(w);
        ASSERT_TRUE(
            near_eq(w, create_vector(5, {blas_norm2(w), 0, 0, 0, 0}), 1e-6));
    }
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
        ASSERT_TRUE(near_eq(
            A.col(0), create_vector(3, {-blas_norm2(A.col(0)), 0, 0}), 1e-6));
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
        auto QR = A;
        vector tau(A.n());
        qr_decomp_hr(QR, tau);
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
    {
        auto A = create_matrix(3, 3,
                               {

                                   12, -51, -12,  //
                                   6, 167, -6,    //
                                   -4, 24, 4,     //
                               });
        auto QR = A;
        vector tau(A.n());
        qr_decomp_hr(QR, tau);
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
        qr_decomp_hr(QR, tau);
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
}
