/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file InvTxfm2dAsmTest.c
 *
 * @brief Unit test for SAD functions:
 * - nxm_sad_kernel_sub_sampled_func
 * - nxm_sad_kernel_func
 * - nxm_sad_averaging_kernel_func
 *
 * @author Cidana-Ryan
 *
 ******************************************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <new>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbComputeSAD.h"
#include "random.h"
#include "util.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random
namespace {
#define MAX_BLOCK_SIZE (128 * 256)
typedef std::tuple<int, int> SadSize;
typedef enum { REF_MAX, SRC_MAX, RANDOM, UNALIGN } TestPattern;

SadSize TEST_SAD_SIZES[] = {SadSize(64, 64),
                            SadSize(64, 32),
                            SadSize(32, 64),
                            SadSize(32, 32),
                            SadSize(32, 16),
                            SadSize(16, 32),
                            SadSize(16, 16),
                            SadSize(16, 8),
                            SadSize(8, 16),
                            SadSize(8, 8),
                            SadSize(8, 4),
                            SadSize(4, 8),
                            SadSize(4, 16),
                            SadSize(16, 4),
                            SadSize(8, 32),
                            SadSize(32, 8),
                            SadSize(16, 64),
                            SadSize(64, 16)};

SadSize TEST_LOOP_SIZES[] = {
    SadSize(16, 16),
    SadSize(16, 8),
    SadSize(8, 16),
    SadSize(8, 8),
    SadSize(8, 4),
    SadSize(4, 8),
    SadSize(4, 16),
    SadSize(16, 4),
};

TestPattern TEST_PATTERNS[] = {REF_MAX, SRC_MAX, RANDOM, UNALIGN};

typedef std::tuple<TestPattern, SadSize> SadTestParam;
/**
 * @brief Base class for SAD test. SADTestBase handle test vector in memory,
 * provide SAD and SAD avg reference function
 *
 */
class SADTestBase : public ::testing::Test {
  public:
    SADTestBase(const int width, const int height, TestPattern test_pattern) {
        width_ = width;
        height_ = height;
        source_stride_ = ref1_stride_ = ref2_stride_ = width_ * 2;
        test_pattern_ = test_pattern;
    };

    void SetUp() override {
        int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (test_pattern_) {
        case REF_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = 0;
                ref1_data_[i] = ref2_data_[i] = mask;
            }
            break;
        }
        case SRC_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = mask;
                ref1_data_[i] = ref2_data_[i] = 0;
            }
            break;
        }
        case RANDOM: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = rnd.random();
                ref1_data_[i] = rnd.random();
                ref2_data_[i] = rnd.random();
            }
            break;
        };
        case UNALIGN: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = rnd.random();
                ref1_data_[i] = rnd.random();
                ref2_data_[i] = rnd.random();
            }
            ref1_stride_ -= 1;
            ref2_stride_ -= 1;
            break;
        }
        default: break;
        }
    }

    uint32_t reference_sad() {
        unsigned int sad = 0;
        for (int h = 0; h < height_; ++h) {
            for (int w = 0; w < width_; ++w) {
                sad += abs(source_data_[h * source_stride_ + w] -
                           ref1_data_[h * ref1_stride_ + w]);
            }
        }
        return sad;
    }

    unsigned int reference_sad_avg() {
        unsigned int sad = 0;
        for (int h = 0; h < height_; ++h) {
            for (int w = 0; w < width_; ++w) {
                const int tmp = ref2_data_[h * ref2_stride_ + w] +
                                ref1_data_[h * ref1_stride_ + w];
                const uint8_t comp_pred = ROUND_POWER_OF_TWO(tmp, 1);
                sad += abs(source_data_[h * source_stride_ + w] - comp_pred);
            }
        }
        return sad;
    }

  protected:
    int width_, height_;
    int source_stride_;
    int ref1_stride_;
    int ref2_stride_;
    TestPattern test_pattern_;

    DECLARE_ALIGNED(32, uint8_t, source_data_[MAX_BLOCK_SIZE]);
    DECLARE_ALIGNED(32, uint8_t, ref1_data_[MAX_BLOCK_SIZE]);
    DECLARE_ALIGNED(32, uint8_t, ref2_data_[MAX_BLOCK_SIZE]);
};
/**
 * @brief Unit test for SAD sub smaple functions include:
 *  - nxm_sad_kernel_sub_sampled_func_ptr_array
 *  - nxm_sad_kernel_sub_sampled_func_ptr_array
 *
 * Test strategy:
 *  This test case combine different wight(4-64) x height(4-64), different test
 * vecotr pattern(MaxRef, MaxSrc, Random, Unalign) to generate test vector. And
 * select correct function from funcPtrArray, run func with test vector, compare
 * result between reference function, non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results come from reference functon, non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *  All functions inside nxm_sad_kernel_sub_sampled_func_ptr_array and
 * nxm_sad_kernel_sub_sampled_func_ptr_array.
 *
 * Test cases:
 *  Width (4, 8, 16, 24, 32, 48, 64) x height 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern (MaxRef, MaxSrc, Random, Unalign).
 *
 */
class SADTestSubSample : public ::testing::WithParamInterface<SadTestParam>,
                         public SADTestBase {
  public:
    SADTestSubSample()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_sad() {
        EbSadKernelNxMType non_avx2_func =
            nxm_sad_kernel_sub_sampled_func_ptr_array[ASM_NON_AVX2]
                                                     [width_ >> 3];
        EbSadKernelNxMType avx2_func =
            nxm_sad_kernel_sub_sampled_func_ptr_array[ASM_AVX2][width_ >> 3];
        ;

        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;
        ref_sad = reference_sad();
        if (non_avx2_func != nullptr)
            non_avx2_sad = non_avx2_func(source_data_,
                                         source_stride_,
                                         ref1_data_,
                                         ref1_stride_,
                                         height_,
                                         width_);
        if (avx2_func != nullptr)
            avx2_sad = avx2_func(source_data_,
                                 source_stride_,
                                 ref1_data_,
                                 ref1_stride_,
                                 height_,
                                 width_);

        if (non_avx2_func != nullptr)
            EXPECT_EQ(ref_sad, non_avx2_sad)
                << "compare ref and non_avx2 error";
        if (avx2_func != nullptr)
            EXPECT_EQ(ref_sad, avx2_sad) << "compare ref and avx2 error";
        if (non_avx2_func != nullptr && avx2_func != nullptr)
            EXPECT_EQ(non_avx2_sad, avx2_sad)
                << "compare non_avx2 and non_avx2 error";
    }
};

TEST_P(SADTestSubSample, SADTestSubSample) {
    check_sad();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADTestSubSample,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_SIZES)));
/**
 * @brief Unit test for SAD functions include:
 *  - nxm_sad_kernel_func_ptr_array
 *  - nxm_sad_kernel_func_ptr_array
 *
 * Test strategy:
 *  This test case combine different wight(4-64) x height(4-64), different test
 * vecotr pattern(MaxRef, MaxSrc, Random, Unalign) to generate test vector. And
 * select correct function from funcPtrArray, run func with test vector, compare
 * result between reference function, non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results come from reference functon, non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *  All functions inside nxm_sad_kernel_func_ptr_array and
 * nxm_sad_kernel_func_ptr_array.
 *
 * Test cases:
 *  Width (4, 8, 16, 24, 32, 48, 64) x height 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern (MaxRef, MaxSrc, Random, Unalign).
 *
 */
class SADTest : public ::testing::WithParamInterface<SadTestParam>,
                public SADTestBase {
  public:
    SADTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_sad() {
        EbSadKernelNxMType non_avx2_func =
            nxm_sad_kernel_func_ptr_array[ASM_NON_AVX2][width_ >> 3];
        EbSadKernelNxMType avx2_func =
            nxm_sad_kernel_func_ptr_array[ASM_AVX2][width_ >> 3];

        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;
        ref_sad = reference_sad();
        if (non_avx2_func != nullptr)
            non_avx2_sad = non_avx2_func(source_data_,
                                         source_stride_,
                                         ref1_data_,
                                         ref1_stride_,
                                         height_,
                                         width_);
        if (avx2_func != nullptr)
            avx2_sad = avx2_func(source_data_,
                                 source_stride_,
                                 ref1_data_,
                                 ref1_stride_,
                                 height_,
                                 width_);

        if (non_avx2_func != nullptr)
            EXPECT_EQ(ref_sad, non_avx2_sad)
                << "compare ref and non_avx2 error";
        if (avx2_func != nullptr)
            EXPECT_EQ(ref_sad, avx2_sad) << "compare ref and avx2 error";
        if (non_avx2_func != nullptr && avx2_func != nullptr)
            EXPECT_EQ(non_avx2_sad, avx2_sad)
                << "compare non_avx2 and non_avx2 error";
    }
};

TEST_P(SADTest, SADTest) {
    check_sad();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_SIZES)));
/**
 * @brief Unit test for SAD Avg functions include:
 *  - nxm_sad_averaging_kernel_func_ptr_array
 *  - nxm_sad_averaging_kernel_func_ptr_array
 *
 * Test strategy:
 *  This test case combine different wight(4-64) x height(4-64), different test
 * vecotr pattern(MaxRef, MaxSrc, Random, Unalign) to generate test vector. And
 * select correct function from funcPtrArray, run func with test vector, compare
 * result between reference function, non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results come from reference functon, non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *  All functions inside nxm_sad_averaging_kernel_func_ptr_array and
 * nxm_sad_averaging_kernel_func_ptr_array.
 *
 * Test cases:
 *  Width (4, 8, 16, 24, 32, 48, 64) x height 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern (MaxRef, MaxSrc, Random, Unalign).
 *
 */
class SADAvgTest : public ::testing::WithParamInterface<SadTestParam>,
                   public SADTestBase {
  public:
    SADAvgTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_sad_avg() {
        EbSadAvgKernelNxMType non_avx2_func =
            nxm_sad_averaging_kernel_func_ptr_array[ASM_NON_AVX2][width_ >> 3];
        EbSadAvgKernelNxMType avx2_func =
            nxm_sad_averaging_kernel_func_ptr_array[ASM_AVX2][width_ >> 3];

        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;
        ref_sad = reference_sad_avg();
        if (non_avx2_func != nullptr)
            non_avx2_sad = non_avx2_func(source_data_,
                                         source_stride_,
                                         ref1_data_,
                                         ref1_stride_,
                                         ref2_data_,
                                         ref2_stride_,
                                         height_,
                                         width_);
        if (avx2_func != nullptr)
            avx2_sad = avx2_func(source_data_,
                                 source_stride_,
                                 ref1_data_,
                                 ref1_stride_,
                                 ref2_data_,
                                 ref2_stride_,
                                 height_,
                                 width_);

        if (non_avx2_func != nullptr)
            EXPECT_EQ(ref_sad, non_avx2_sad)
                << "compare ref and non_avx2 error";
        if (avx2_func != nullptr)
            EXPECT_EQ(ref_sad, avx2_sad) << "compare ref and avx2 error";
        if (non_avx2_func != nullptr && avx2_func != nullptr)
            EXPECT_EQ(non_avx2_sad, avx2_sad)
                << "compare non_avx2 and non_avx2 error";
    }
};

TEST_P(SADAvgTest, SADAvgTest) {
    check_sad_avg();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADAvgTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_SIZES)));
/**
 * @brief Unit test for SAD loop (sparse) functions include:
 *  - sad_loop_kernel_sparse_{sse4_1,avx2}_intrin
 *  - sad_loop_kernel_{sse4_1,avx2}_intrin
 *
 * Test strategy:
 *  This test case combine different wight(4-16) x height(4-16), different test
 * vecotr pattern(MaxRef, MaxSrc, Random, Unalign) to generate test vector.
 * Run func with test vector, compare result between  non_avx2 function and avx2
 * function.
 *
 *
 * Expect result:
 *  Results come from  non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 *
 */
class SadLoopTest : public ::testing::WithParamInterface<SadTestParam>,
                    public SADTestBase {
  public:
    SadLoopTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
        source_stride_ = ref_1_stride_ = ref2_stride_ = 256;
    };

  protected:
    void check_sad_loop_sparse() {
        uint64_t best_sad1[4];
        int16_t x_search_center1 = 0;
        int16_t y_search_center1 = 0;

        printf("Width:%d, height_:%d\r\n", width_, height_);
        sad_loop_kernel_sparse_sse4_1_intrin(source_data_,
                                             source_stride_,
                                             ref_1_data_,
                                             ref_1_stride_,
                                             height_,
                                             width_,
                                             best_sad1,
                                             &x_search_center1,
                                             &y_search_center1,
                                             source_stride_,
                                             128,
                                             128);
        uint64_t best_sad2[4];
        int16_t x_search_center2 = 0;
        int16_t y_search_center2 = 0;
        sad_loop_kernel_sparse_avx2_intrin(source_data_,
                                           source_stride_,
                                           ref_1_data_,
                                           ref_1_stride_,
                                           height_,
                                           width_,
                                           best_sad2,
                                           &x_search_center2,
                                           &y_search_center2,
                                           source_stride_,
                                           128,
                                           128);

        EXPECT_EQ(0, memcmp(best_sad1, best_sad2, sizeof(best_sad1)))
            << "compare bast_sad error";
        EXPECT_EQ(x_search_center1, x_search_center2)
            << "compare x_search_center error";
        EXPECT_EQ(y_search_center1, y_search_center2)
            << "compare y_search_center error";
    }

    void check_sad_loop() {
        uint64_t best_sad1[4];
        int16_t x_search_center1 = 0;
        int16_t y_search_center1 = 0;

        sad_loop_kernel_sse4_1_intrin(source_data_,
                                      source_stride_,
                                      ref_1_data_,
                                      ref_1_stride_,
                                      16,
                                      16,
                                      best_sad1,
                                      &x_search_center1,
                                      &y_search_center1,
                                      source_stride_,
                                      128,
                                      128);
        uint64_t best_sad2[4];
        int16_t x_search_center2 = 0;
        int16_t y_search_center2 = 0;
        sad_loop_kernel_avx2_intrin(source_data_,
                                    source_stride_,
                                    ref_1_data_,
                                    ref_1_stride_,
                                    16,
                                    16,
                                    best_sad2,
                                    &x_search_center2,
                                    &y_search_center2,
                                    source_stride_,
                                    128,
                                    128);

        EXPECT_EQ(0, memcmp(best_sad1, best_sad2, sizeof(best_sad1)))
            << "compare bast_sad error";
        EXPECT_EQ(x_search_center1, x_search_center2)
            << "compare x_search_center error";
        EXPECT_EQ(y_search_center1, y_search_center2)
            << "compare y_search_center error";
    }

    int source_stride_;
    int ref_1_stride_;
    TestPattern vector_param_;

    DECLARE_ALIGNED(32, uint8_t, source_data_[MAX_BLOCK_SIZE]);
    DECLARE_ALIGNED(32, uint8_t, ref_1_data_[MAX_BLOCK_SIZE]);
};

TEST_P(SadLoopTest, SadLoopSparseTest) {
    check_sad_loop_sparse();
}

TEST_P(SadLoopTest, SadLoopTest) {
    check_sad_loop();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SadLoopTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_LOOP_SIZES)));

/**
 * @brief Unit test for GetEightSadTest functions include:
 *  -
 * sget_eight_horizontal_search_point_results_8x8_16x16_pu_{sse41,avx2}_intrin
 *
 * Test strategy:
 *  This test use different test vecotr pattern(MaxRef, MaxSrc, Random, Unalign)
 * to generate test vector. Run func with test vector, compare result between
 * non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results come from  non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 *
 */

class GetEightSadTest : public ::testing::TestWithParam<TestPattern> {
  public:
    GetEightSadTest() : vector_param_(GetParam()) {
        source_stride_ = ref_1_stride_ = 256;
    };
    void SetUp() override {
        int32_t mask = (1 << 8) - 1;
        SVTRandom *rnd = new SVTRandom(0, mask);
        switch (vector_param_) {
        case REF_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = 0;
                ref_1_data_[i] = mask;
            }
            break;
        }
        case SRC_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = mask;
                ref_1_data_[i] = 0;
            }
            break;
        }
        case RANDOM: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = rnd->random();
                ref_1_data_[i] = rnd->random();
            }
        } break;
        case UNALIGN: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                source_data_[i] = rnd->random();
                ref_1_data_[i] = rnd->random();
            }
            source_stride_ -= 1;
            ref_1_stride_ -= 1;
        } break;
        default: break;
        }
        delete rnd;
    }

  protected:
    void check_get_eight() {
        uint32_t best_sad8x8_1[4], best_mv8x8_1[4], best_sad16x16_1,
            best_mv16x16_1;
        uint16_t sad16x16_1[8];
        get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin(
            source_data_,
            source_stride_,
            ref_1_data_,
            ref_1_stride_,
            best_sad8x8_1,
            best_mv8x8_1,
            &best_sad16x16_1,
            &best_mv16x16_1,
            0,
            sad16x16_1,
            false);

        uint32_t best_sad8x8_2[4], best_mv8x8_2[4], best_sad16x16_2,
            best_mv16x16_2;
        uint16_t sad16x16_2[8];
        get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin(
            source_data_,
            source_stride_,
            ref_1_data_,
            ref_1_stride_,
            best_sad8x8_2,
            best_mv8x8_2,
            &best_sad16x16_2,
            &best_mv16x16_2,
            0,
            sad16x16_2,
            false);
        EXPECT_EQ(0, memcmp(best_sad8x8_1, best_sad8x8_2, 4));
        EXPECT_EQ(0, memcmp(best_mv8x8_1, best_mv8x8_2, 4));
        EXPECT_EQ(best_sad16x16_1, best_sad16x16_2);
        EXPECT_EQ(best_mv16x16_1, best_mv16x16_2);
        EXPECT_EQ(sad16x16_1, sad16x16_2);
    }
    int source_stride_;
    int ref_1_stride_;
    TestPattern vector_param_;

    DECLARE_ALIGNED(32, uint8_t, source_data_[MAX_BLOCK_SIZE]);
    DECLARE_ALIGNED(32, uint8_t, ref_1_data_[MAX_BLOCK_SIZE]);
};

TEST_P(GetEightSadTest, GetEightSadTest) {
    check_get_eight();
}

INSTANTIATE_TEST_CASE_P(SAD, GetEightSadTest,
                        ::testing::ValuesIn(TEST_PATTERNS));

}  // namespace
