#include <iostream>
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "matrix.hpp"
#include "matrix_oper.hpp"

using namespace markovgg;

vec create_vec(const std::vector<real_t>& val)
{
    vec v(val.size(), 0.0);
    for (size_t i = 0; i < v.dim(); i++)
    {
        v[i] = val[i];
    }
    return v;
}

sqr_mat create_sqr_mat(size_t dim, const std::vector<real_t>& val)
{
    sqr_mat_crtor mat_crt(dim);
    assert(dim * dim == val.size());
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            mat_crt.add_entry(i, j, val[i * dim + j]);
        }
    }
    return mat_crt.create();
}

matrix create_matrix(size_t row_dim, size_t col_dim,
                     const std::vector<real_t>& val)
{
    matrix_crtor mat_crt(row_dim, col_dim);
    assert(row_dim * col_dim == val.size());
    for (size_t i = 0; i < row_dim; i++)
    {
        for (size_t j = 0; j < col_dim; j++)
        {
            mat_crt.add_entry(i, j, val[i * col_dim + j]);
        }
    }
    return mat_crt.create();
}

TEST(test_matrix, sqr_create)
{
    sqr_mat_crtor mat_crt(3);
    mat_crt.add_entry(2, 2, 3);
    mat_crt.add_entry(0, 0, 3);
    mat_crt.add_entry(1, 2, 1.5);
    mat_crt.add_entry(1, 1, 1);
    mat_crt.add_entry(1, 1, 1);
    mat_crt.add_entry(1, 2, 1.5);
    mat_crt.add_entry(1, 1, 1);
    sqr_mat cmat = mat_crt.create();
    print_mat(cmat);
    ASSERT_EQ(cmat(1, 2), 3.0);
    ASSERT_EQ(cmat(2, 2), 3.0);
    ASSERT_EQ(cmat(1, 0), 0.0);
}

TEST(test_matrix, matrix_create)
{
    matrix_crtor mat_crt(3, 4);
    mat_crt.add_entry(0, 0, 1.5);
    mat_crt.add_entry(0, 0, 1.5);
    mat_crt.add_entry(0, 1, 2.5);
    mat_crt.add_entry(0, 1, 2.5);
    mat_crt.add_entry(1, 1, 3.5);
    mat_crt.add_entry(1, 1, 3.5);
    mat_crt.add_entry(1, 2, 4.5);
    mat_crt.add_entry(1, 2, 4.5);
    mat_crt.add_entry(2, 2, 5.5);
    mat_crt.add_entry(2, 2, 5.5);
    mat_crt.add_entry(2, 3, 6.5);
    mat_crt.add_entry(2, 3, 6.5);
    auto mat = mat_crt.create();
    print_mat(mat);
    ASSERT_EQ(mat(0, 0), 3);
    ASSERT_EQ(mat(0, 1), 5);
    ASSERT_EQ(mat(0, 2), 0);
    ASSERT_EQ(mat(0, 3), 0);
    ASSERT_EQ(mat(1, 0), 0);
    ASSERT_EQ(mat(1, 1), 7);
    ASSERT_EQ(mat(1, 2), 9);
    ASSERT_EQ(mat(1, 3), 0);
    ASSERT_EQ(mat(2, 0), 0);
    ASSERT_EQ(mat(2, 1), 0);
    ASSERT_EQ(mat(2, 2), 11);
    ASSERT_EQ(mat(2, 3), 13);
}

TEST(test_matrix, partition_sqr)
{
    auto mat = create_sqr_mat(5, {1,  2,  3,  4,  5,

                                  6,  7,  8,  9,  10,

                                  11, 12, 13, 14, 15,

                                  16, 17, 18, 19, 20,

                                  21, 22, 23, 24, 25});
    print_mat(mat);
    auto part0 = part_mat(mat, {1, 3, 4});
    print_mat(part0);
    ASSERT_EQ(part0(0, 0), 7);
    ASSERT_EQ(part0(2, 2), 25);
    auto part1 = part_mat(mat, {1, 3, 4}, {2, 3});
    print_mat(part1);
    ASSERT_EQ(part1(0, 0), 8);
    ASSERT_EQ(part1(0, 1), 9);
    ASSERT_EQ(part1(2, 1), 24);
}

TEST(test_matrix, partition_matrix)
{
    auto mat = create_matrix(5, 5, {1,  2,  3,  4,  5,

                                    6,  7,  8,  9,  10,

                                    11, 12, 13, 14, 15,

                                    16, 17, 18, 19, 20,

                                    21, 22, 23, 24, 25});
    auto part0 = part_mat(mat, {1, 3, 4});
    print_mat(part0);
    ASSERT_EQ(part0(0, 0), 7);
    ASSERT_EQ(part0(2, 2), 25);
    auto part1 = part_mat(mat, {1, 3, 4}, {2, 3});
    print_mat(part1);
    ASSERT_EQ(part1(0, 0), 8);
    ASSERT_EQ(part1(0, 1), 9);
    ASSERT_EQ(part1(2, 1), 24);
}

TEST(test_matrix, partition_vec)
{
    auto v = create_vec({0, 1, 2, 3, 4, 5});
    auto part_v = part_vec(v, {2, 3, 5});
    print_vec(part_v);
    ASSERT_TRUE(near_eq(part_v, create_vec({2, 3, 5}), 1e-20));
}

TEST(test_matrix, mul)
{
    print_sep();
    auto mat = create_matrix(4, 3,
                             {1, 2, 3,

                              4, 5, 6,

                              7, 8, 9,

                              10, 11, 12});
    auto left = create_vec({1, 2, 3, 4});
    auto ans = create_vec({70, 80, 90});
    print_vec(left);
    std::cout << "* ";
    print_mat(mat);
    std::cout << "= ";
    auto right = mul(left, mat);
    print_vec(right);
    ASSERT_EQ(ans, right);
    acc_mul(right, left, mat);
    print_vec(right);
    ASSERT_EQ(mul(ans, 2), right);
}
