#include "blas.hpp"
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "linalg.hpp"
#include "matvec_oper.hpp"
#include "splinalg.hpp"
#include "spmatrix_oper.hpp"

static real_t tol = 1e-6;

TEST(test_large, rgmres_example)
{
    std::ifstream matf("./test/matrix/example/example.mtx");
    std::ifstream bf("./test/matrix/example/example_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    vector x(b.dim(), 0.0);
    real_t prec = spsolve_restart_gmres_gms(A, x, b, 30, tol, 1000, 30);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}

TEST(test_large, rgmres_circuit_1)
{
    std::ifstream matf("./test/matrix/circuit_1/circuit_1.mtx");
    std::ifstream bf("./test/matrix/circuit_1/circuit_1_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    vector x(b.dim(), 0.0);
    real_t prec = spsolve_restart_gmres_gms(A, x, b, 30, tol, 100, 30);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}

TEST(test_large, rgmres_chem_master1)
{
    std::ifstream matf("./test/matrix/chem_master1/chem_master1.mtx");
    std::ifstream bf("./test/matrix/chem_master1/chem_master1_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    vector x(b.dim(), 0.0);
    real_t prec = spsolve_restart_gmres_gms(A, x, b, 100, tol, 100, 10);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}

TEST(test_large, sor_chem_master1)
{
    std::ifstream matf("./test/matrix/chem_master1/chem_master1.mtx");
    std::ifstream bf("./test/matrix/chem_master1/chem_master1_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    auto x = vector(b.dim(), 0.0);
    real_t w = 1.2;
    real_t prec = spsolve_sor_method(A, x, b, w, tol, 10000, 100);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}

TEST(test_large, lu_circuit_1)
{
    std::ifstream matf("./test/matrix/circuit_1/circuit_1.mtx");
    std::ifstream bf("./test/matrix/circuit_1/circuit_1_b.mtx");
    const matrix A_ = spmatrix2dense(read_cord_format(matf));
    const vector b_ = vector(read_array_format(bf).col(0));
    matrix U = A_;
    auto P = identity_matrix(A_.m());
    decomp_lup(U, P);
    matrix L(A_.m(), A_.n(), 0.0);
    unpack_lu(U, L);
    ASSERT_TRUE(near_eq(L * U, P * A_, 1e-6));
}

TEST(test_large, solve_lu_circuit_1)
{
    std::ifstream matf("./test/matrix/circuit_1/circuit_1.mtx");
    std::ifstream bf("./test/matrix/circuit_1/circuit_1_b.mtx");
    const matrix A_ = spmatrix2dense(read_cord_format(matf));
    const vector b_ = vector(read_array_format(bf).col(0));
    matrix A = A_;
    vector x = b_;
    solve_lu(A, x);
    ASSERT_TRUE(near_eq(A_ * x, b_, tol));
}

TEST(test_large, solve_qr_circuit_1)
{
    std::ifstream matf("./test/matrix/circuit_1/circuit_1.mtx");
    std::ifstream bf("./test/matrix/circuit_1/circuit_1_b.mtx");
    const matrix A_ = spmatrix2dense(read_cord_format(matf));
    const vector b_ = vector(read_array_format(bf).col(0));
    matrix A = A_;
    vector x = b_;
    solve_qr(A, x);
    ASSERT_TRUE(near_eq(A_ * x, b_, tol));
}

TEST(test_large, precond_rgmres_chem_master1)
{
    std::ifstream matf("./test/matrix/chem_master1/chem_master1.mtx");
    std::ifstream bf("./test/matrix/chem_master1/chem_master1_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    spmatrix iLU = A;
    spdecomp_ilu(iLU);
    vector x(b.dim(), 0.0);
    real_t prec = spsolve_restart_gmres_gms(
        A, x, b, 30, tol, 100, 10,
        [&](vector_mutable_view x) -> void { spsolve_ilu(iLU, x); });

    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}
