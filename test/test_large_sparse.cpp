#include "blas.hpp"
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "matvec_oper.hpp"
#include "splinalg.hpp"
#include "spmatrix_oper.hpp"

static real_t tol = 1e-10;

TEST(test_large_sparse, rgmres_example)
{
    std::ifstream matf("./test/matrix/example/example.mtx");
    std::ifstream bf("./test/matrix/example/example_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    vector x(b.dim(), 0.0);
    real_t prec = linsv_restart_gmres_gms(A, x, b, 30, tol, 1000, 30);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}

TEST(test_large_sparse, rgmres_circuit_1)
{
    std::ifstream matf("./test/matrix/circuit_1/circuit_1.mtx");
    std::ifstream bf("./test/matrix/circuit_1/circuit_1_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    vector x(b.dim(), 0.0);
    real_t prec = linsv_restart_gmres_gms(A, x, b, 30, tol, 100, 30);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}

TEST(test_large_sparse, sor_chem_master1)
{
    std::ifstream matf("./test/matrix/chem_master1/chem_master1.mtx");
    std::ifstream bf("./test/matrix/chem_master1/chem_master1_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    auto x = vector(b.dim(), 0.0);
    real_t w = 1.2;
    linsv_sor_method(A, x, b, w, tol, 1000, 100);
    ASSERT_TRUE(near_eq(dot(A, x), b, 1e-4));
}

TEST(test_large_sparse, rgmres_chem_master1)
{
    std::ifstream matf("./test/matrix/chem_master1/chem_master1.mtx");
    std::ifstream bf("./test/matrix/chem_master1/chem_master1_b.mtx");
    const spmatrix A = read_cord_format(matf);
    const vector b = vector(read_array_format(bf).col(0));
    vector x(b.dim(), 0.0);
    real_t prec = linsv_restart_gmres_gms(A, x, b, 100, tol, 100, 20);
    ASSERT_LT(prec, tol);
    ASSERT_TRUE(near_eq(dot(A, x), b, tol));
}
