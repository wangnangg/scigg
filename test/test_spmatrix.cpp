#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "spmatrix.hpp"
#include "spmatvec_oper.hpp"

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

TEST(test_spmatrix, spmatrix2dense_crs)
{
    spmatrix spm = create_spmatrix(5, 4,
                                   {
                                       1,  2,  3,  4,   //
                                       5,  6,  7,  8,   //
                                       9,  10, 11, 12,  //
                                       13, 14, 15, 16,  //
                                       17, 18, 19, 20   //
                                   },
                                   true);
    print(spm);
    matrix dsm = spmatrix2dense(spm);
    print(dsm);
    ASSERT_EQ(dsm, create_matrix(5, 4,
                                 {
                                     1,  2,  3,  4,   //
                                     5,  6,  7,  8,   //
                                     9,  10, 11, 12,  //
                                     13, 14, 15, 16,  //
                                     17, 18, 19, 20   //
                                 }

                                 ));
}

TEST(test_spmatrix, spmatrix2dense_ccs)
{
    spmatrix spm = create_spmatrix(5, 4,
                                   {
                                       1,  2,  3,  4,   //
                                       5,  6,  7,  8,   //
                                       9,  10, 11, 12,  //
                                       13, 14, 15, 16,  //
                                       17, 18, 19, 20   //
                                   },
                                   false);
    print(spm);
    matrix dsm = spmatrix2dense(spm);
    print(dsm);
    ASSERT_EQ(dsm, create_matrix(5, 4,
                                 {
                                     1,  2,  3,  4,   //
                                     5,  6,  7,  8,   //
                                     9,  10, 11, 12,  //
                                     13, 14, 15, 16,  //
                                     17, 18, 19, 20   //
                                 }

                                 ));
}
