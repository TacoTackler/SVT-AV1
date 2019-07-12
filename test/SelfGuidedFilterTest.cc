/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file AdaptiveScanTest.cc
 *
 * @brief Unit test for self guided filter
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/

#include "gtest/gtest.h"
#include <stdlib.h>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbDefinitions.h"
#include "random.h"
#include "EbRestoration.h"
#include "EbUtility.h"
#include "aom_dsp_rtcd.h"
#include "util.h"

namespace {
using ::testing::make_tuple;
using ::testing::tuple;
using svt_av1_test_tool::SVTRandom;

using SgrFunc = void (*)(const uint8_t *dat8, int width, int height, int stride,
                         int eps, const int *xqd, uint8_t *dst8, int dst_stride,
                         int32_t *tmpbuf, int bit_depth, int highbd);

typedef tuple<SgrFunc, int> FilterTestParam;
template <typename Sample>
class SelfGuidedFilterTest : public ::testing::TestWithParam<FilterTestParam> {
  public:
    virtual ~SelfGuidedFilterTest() {
    }

    virtual void SetUp() {
        input_data_ = (Sample *)aom_memalign(
            32, in_stride_ * (max_h_ + 32) * sizeof(Sample));
        output_tst_data_ = (Sample *)aom_memalign(
            32, out_stride_ * (max_h_ + 32) * sizeof(Sample));
        output_ref_data_ = (Sample *)aom_memalign(
            32, out_stride_ * (max_h_ + 32) * sizeof(Sample));
        tmpbuf_ = (int32_t *)aom_memalign(32, RESTORATION_TMPBUF_SIZE);

        input_ = input_data_ + in_stride_ * 16 + 16;
        output_tst_ = output_tst_data_ + out_stride_ * 16 + 16;
        output_ref_ = output_ref_data_ + out_stride_ * 16 + 16;
        ref_fun_ = apply_selfguided_restoration_c;
    }

    virtual void TearDown() {
        aom_free(input_data_);
        aom_free(output_tst_data_);
        aom_free(output_ref_data_);
        aom_free(tmpbuf_);
        aom_clear_system_state();
    }

  protected:
    void RunCorrectnessTest() {
        tst_fun_ = TEST_GET_PARAM(0);
        bit_depth_ = TEST_GET_PARAM(1);

        int i, j, k;
        SVTRandom rnd(bit_depth_, false);
        for (i = 0; i < NUM_ITERS_; ++i) {
            // prepare random data
            for (j = -16; j < max_h_ + 16; ++j)
                for (k = -16; k < max_w_ + 16; ++k)
                    input_[j * in_stride_ + k] = rnd.random();

            SVTRandom xqd_rnd1(0, SGRPROJ_PRJ_MAX0 + 1 - SGRPROJ_PRJ_MIN0);
            SVTRandom xqd_rnd2(0, SGRPROJ_PRJ_MAX1 + 1 - SGRPROJ_PRJ_MIN1);
            SVTRandom eps_rnd(0, 1 << SGRPROJ_PARAMS_BITS);
            int xqd[2] = {SGRPROJ_PRJ_MIN0 + xqd_rnd1.random(),
                          SGRPROJ_PRJ_MIN1 + xqd_rnd2.random()};
            int eps = eps_rnd.random();

            // Test various tile sizes around 256x256
            int test_w = max_w_ - (i / 9);
            int test_h = max_h_ - (i % 9);

            for (k = 0; k < test_h; k += pu_height_) {
                for (j = 0; j < test_w; j += pu_width_) {
                    int w = AOMMIN(pu_width_, test_w - j);
                    int h = AOMMIN(pu_height_, test_h - k);
                    run_filter(j, k, w, h, eps, xqd);
                }
            }

            for (j = 0; j < test_h; ++j) {
                for (k = 0; k < test_w; ++k) {
                    ASSERT_EQ(output_tst_[j * out_stride_ + k],
                              output_ref_[j * out_stride_ + k]);
                }
            }
        }
    }

    virtual void run_filter(int offset_x, int offset_y, int width, int height,
                            int eps, int *xqd) = 0;

  protected:
    SgrFunc tst_fun_;
    SgrFunc ref_fun_;
    int bit_depth_;
    Sample *input_data_;
    Sample *output_tst_data_;
    Sample *output_ref_data_;
    int32_t *tmpbuf_;
    Sample *input_;
    Sample *output_tst_;
    Sample *output_ref_;
    static const int pu_width_ = RESTORATION_PROC_UNIT_SIZE;
    static const int pu_height_ = RESTORATION_PROC_UNIT_SIZE;
    // Set the maximum width/height to test here. We actually test a small
    // range of sizes *up to* this size, so that we can check, eg.,
    // the behaviour on tiles which are not a multiple of 4 wide.
    static const int max_w_ = 260, max_h_ = 260, in_stride_ = 672,
                     out_stride_ = 672;
    static const int NUM_ITERS_ = 81;
};

class LbdSelfGuidedFilterTest : public SelfGuidedFilterTest<uint8_t> {
    void run_filter(int offset_x, int offset_y, int width, int height, int eps,
                    int *xqd) override {
        uint8_t *input_p = input_ + offset_y * in_stride_ + offset_x;
        uint8_t *output_tst_p = output_tst_ + offset_y * out_stride_ + offset_x;
        uint8_t *output_ref_p = output_ref_ + offset_y * out_stride_ + offset_x;
        tst_fun_(input_p,
                 width,
                 height,
                 in_stride_,
                 eps,
                 xqd,
                 output_tst_p,
                 out_stride_,
                 tmpbuf_,
                 bit_depth_,
                 0);
        ref_fun_(input_p,
                 width,
                 height,
                 in_stride_,
                 eps,
                 xqd,
                 output_ref_p,
                 out_stride_,
                 tmpbuf_,
                 bit_depth_,
                 0);
    }
};

TEST_P(LbdSelfGuidedFilterTest, CorrectnessTest) {
    RunCorrectnessTest();
}

INSTANTIATE_TEST_CASE_P(
    DeblockTest, LbdSelfGuidedFilterTest,
    ::testing::Combine(::testing::Values(apply_selfguided_restoration_avx2),
                       ::testing::Values(8)));

class HbdSelfGuidedFilterTest : public SelfGuidedFilterTest<uint16_t> {
    void run_filter(int offset_x, int offset_y, int width, int height, int eps,
                    int *xqd) override {
        uint16_t *input_p = input_ + offset_y * in_stride_ + offset_x;
        uint16_t *output_tst_p =
            output_tst_ + offset_y * out_stride_ + offset_x;
        uint16_t *output_ref_p =
            output_ref_ + offset_y * out_stride_ + offset_x;
        tst_fun_(CONVERT_TO_BYTEPTR(input_p),
                 width,
                 height,
                 in_stride_,
                 eps,
                 xqd,
                 CONVERT_TO_BYTEPTR(output_tst_p),
                 out_stride_,
                 tmpbuf_,
                 bit_depth_,
                 1);
        ref_fun_(CONVERT_TO_BYTEPTR(input_p),
                 width,
                 height,
                 in_stride_,
                 eps,
                 xqd,
                 CONVERT_TO_BYTEPTR(output_ref_p),
                 out_stride_,
                 tmpbuf_,
                 bit_depth_,
                 1);
    }
};

TEST_P(HbdSelfGuidedFilterTest, CorrectnessTest) {
    RunCorrectnessTest();
}

const int highbd_params_avx2[] = {8, 10, 12};
INSTANTIATE_TEST_CASE_P(
    DeblockTest, HbdSelfGuidedFilterTest,
    ::testing::Combine(::testing::Values(apply_selfguided_restoration_avx2),
                       ::testing::ValuesIn(highbd_params_avx2)));
}  // namespace
