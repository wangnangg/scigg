#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "spmatrix.hpp"

TEST(test_spmatrix, create)
{
    spmatrix_creator crtor(3, 5);
    crtor.add_entry(2, 3, 1.0);
    crtor.add_entry(1, 2, 2.0);
    auto mat1 = crtor.create(true);
    auto mat2 = crtor.create(false);
    print(mat1);
    print(mat2);
    ASSERT_EQ(mat1(2, 3), 1.0);
    ASSERT_EQ(mat1(1, 2), 2.0);
    ASSERT_EQ(mat1(0, 0), 0.0);
    ASSERT_EQ(mat2(2, 3), 1.0);
    ASSERT_EQ(mat2(1, 2), 2.0);
    ASSERT_EQ(mat2(0, 0), 0.0);
}
