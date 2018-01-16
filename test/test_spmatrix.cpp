#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
#include "matvec_oper.hpp"
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

TEST(test_spmatrix, spmatrix_add1)
{
    spmatrix spm1 = create_spmatrix(5, 4,
                                    {
                                        1,  2,  3,  4,   //
                                        5,  6,  7,  8,   //
                                        9,  10, 11, 12,  //
                                        13, 14, 15, 16,  //
                                        17, 18, 19, 20   //
                                    },
                                    false);
    spmatrix spm2 = create_spmatrix(5, 4,
                                    {
                                        1 * 2,  2 * 2,  3 * 2,  4 * 2,   //
                                        5 * 2,  6 * 2,  7 * 2,  8 * 2,   //
                                        9 * 2,  10 * 2, 11 * 2, 12 * 2,  //
                                        13 * 2, 14 * 2, 15 * 2, 16 * 2,  //
                                        17 * 2, 18 * 2, 19 * 2, 20 * 2   //
                                    },
                                    false);
    spmatrix spm3 = add(spm1, spm2);
    ASSERT_EQ(spmatrix2dense(spm3),
              create_matrix(5, 4,
                            {
                                1 * 3,  2 * 3,  3 * 3,  4 * 3,   //
                                5 * 3,  6 * 3,  7 * 3,  8 * 3,   //
                                9 * 3,  10 * 3, 11 * 3, 12 * 3,  //
                                13 * 3, 14 * 3, 15 * 3, 16 * 3,  //
                                17 * 3, 18 * 3, 19 * 3, 20 * 3   //
                            }));
}

TEST(test_spmatrix, spmatrix_add2)
{
    spmatrix spm1 = create_spmatrix(3, 3,
                                    {
                                        1, 0, 0,  //
                                        1, 2, 0,  //
                                        0, 0, 0   //
                                    },
                                    false);
    spmatrix spm2 = create_spmatrix(3, 3,
                                    {
                                        1, 0, 1,  //
                                        0, 2, 1,  //
                                        0, 1, 0   //
                                    },
                                    false);

    ASSERT_EQ(spmatrix2dense(spm1 + spm2),
              spmatrix2dense(spm1) + spmatrix2dense(spm2));
}
