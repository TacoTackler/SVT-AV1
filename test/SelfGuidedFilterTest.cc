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

void av1_loop_restoration_precal() {
#if 0
  GenSgrprojVtable();
#endif
}

typedef void (*SgrFunc)(const uint8_t *dat8, int width, int height, int stride,
                        int eps, const int *xqd, uint8_t *dst8, int dst_stride,
                        int32_t *tmpbuf, int bit_depth, int highbd);

// Test parameter list:
//  <tst_fun_>
typedef tuple<SgrFunc> FilterTestParam;

class AV1SelfguidedFilterTest
    : public ::testing::TestWithParam<FilterTestParam> {
  public:
    virtual ~AV1SelfguidedFilterTest() {
    }
    virtual void SetUp() {
    }

    virtual void TearDown() {
        aom_clear_system_state();
    }

  protected:
    void RunCorrectnessTest() {
        tst_fun_ = TEST_GET_PARAM(0);
        const int pu_width = RESTORATION_PROC_UNIT_SIZE;
        const int pu_height = RESTORATION_PROC_UNIT_SIZE;
        // Set the maximum width/height to test here. We actually test a small
        // range of sizes *up to* this size, so that we can check, eg.,
        // the behaviour on tiles which are not a multiple of 4 wide.
        const int max_w = 260, max_h = 260, stride = 672, out_stride = 672;
        const int NUM_ITERS = 81;
        int i, j, k;
        int bit_depth = 8;

        uint8_t *input_ = (uint8_t *)aom_memalign(
            32, stride * (max_h + 32) * sizeof(uint8_t));
        uint8_t *output_ = (uint8_t *)aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint8_t));
        uint8_t *output2_ = (uint8_t *)aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint8_t));
        int32_t *tmpbuf = (int32_t *)aom_memalign(32, RESTORATION_TMPBUF_SIZE);

        uint8_t *input = input_ + stride * 16 + 16;
        uint8_t *output = output_ + out_stride * 16 + 16;
        uint8_t *output2 = output2_ + out_stride * 16 + 16;

        // ACMRandom rnd(ACMRandom::DeterministicSeed());
        SVTRandom rnd(bit_depth, false);
        av1_loop_restoration_precal();

        for (i = 0; i < NUM_ITERS; ++i) {
            for (j = -16; j < max_h + 16; ++j)
                for (k = -16; k < max_w + 16; ++k)
                    input[j * stride + k] =
                        rnd.random();  // rnd.Rand16() & 0xFF;

            SVTRandom rnd1(0, SGRPROJ_PRJ_MAX0 + 1 - SGRPROJ_PRJ_MIN0);
            SVTRandom rnd2(0, SGRPROJ_PRJ_MAX1 + 1 - SGRPROJ_PRJ_MIN1);
            SVTRandom eps_rnd(0, 1 << SGRPROJ_PARAMS_BITS);
            int xqd[2] = {SGRPROJ_PRJ_MIN0 + rnd1.random(),
                          SGRPROJ_PRJ_MIN1 + rnd2.random()};
            int eps = eps_rnd.random();

            // Test various tile sizes around 256x256
            int test_w = max_w - (i / 9);
            int test_h = max_h - (i % 9);

            for (k = 0; k < test_h; k += pu_height)
                for (j = 0; j < test_w; j += pu_width) {
                    int w = AOMMIN(pu_width, test_w - j);
                    int h = AOMMIN(pu_height, test_h - k);
                    uint8_t *input_p = input + k * stride + j;
                    uint8_t *output_p = output + k * out_stride + j;
                    uint8_t *output2_p = output2 + k * out_stride + j;
                    tst_fun_(input_p,
                             w,
                             h,
                             stride,
                             eps,
                             xqd,
                             output_p,
                             out_stride,
                             tmpbuf,
                             8,
                             0);
                    apply_selfguided_restoration_c(input_p,
                                                   w,
                                                   h,
                                                   stride,
                                                   eps,
                                                   xqd,
                                                   output2_p,
                                                   out_stride,
                                                   tmpbuf,
                                                   8,
                                                   0);
                }

            for (j = 0; j < test_h; ++j)
                for (k = 0; k < test_w; ++k) {
                    ASSERT_EQ(output[j * out_stride + k],
                              output2[j * out_stride + k]);
                }
        }

        aom_free(input_);
        aom_free(output_);
        aom_free(output2_);
        aom_free(tmpbuf);
    }

  private:
    SgrFunc tst_fun_;
};

TEST_P(AV1SelfguidedFilterTest, CorrectnessTest) {
    RunCorrectnessTest();
}

INSTANTIATE_TEST_CASE_P(AVX2, AV1SelfguidedFilterTest,
                        ::testing::Values(apply_selfguided_restoration_avx2));

// Test parameter list:
//  <tst_fun_, bit_depth>
typedef tuple<SgrFunc, int> HighbdFilterTestParam;

class AV1HighbdSelfguidedFilterTest
    : public ::testing::TestWithParam<HighbdFilterTestParam> {
  public:
    virtual ~AV1HighbdSelfguidedFilterTest() {
    }
    virtual void SetUp() {
    }

    virtual void TearDown() {
        aom_clear_system_state();
    }

  protected:
    void RunCorrectnessTest() {
        tst_fun_ = TEST_GET_PARAM(0);
        const int pu_width = RESTORATION_PROC_UNIT_SIZE;
        const int pu_height = RESTORATION_PROC_UNIT_SIZE;
        // Set the maximum width/height to test here. We actually test a small
        // range of sizes *up to* this size, so that we can check, eg.,
        // the behaviour on tiles which are not a multiple of 4 wide.
        const int max_w = 260, max_h = 260, stride = 672, out_stride = 672;
        const int NUM_ITERS = 81;
        int i, j, k;
        int bit_depth = TEST_GET_PARAM(1);

        uint16_t *input_ = (uint16_t *)aom_memalign(
            32, stride * (max_h + 32) * sizeof(uint16_t));
        uint16_t *output_ = (uint16_t *)aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint16_t));
        uint16_t *output2_ = (uint16_t *)aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint16_t));
        int32_t *tmpbuf = (int32_t *)aom_memalign(32, RESTORATION_TMPBUF_SIZE);

        uint16_t *input = input_ + stride * 16 + 16;
        uint16_t *output = output_ + out_stride * 16 + 16;
        uint16_t *output2 = output2_ + out_stride * 16 + 16;

        // ACMRandom rnd(ACMRandom::DeterministicSeed());
        SVTRandom rnd(bit_depth, false);

        av1_loop_restoration_precal();

        for (i = 0; i < NUM_ITERS; ++i) {
            for (j = -16; j < max_h + 16; ++j)
                for (k = -16; k < max_w + 16; ++k)
                    input[j * stride + k] = rnd.random();

#if 0
      int xqd[2] = { SGRPROJ_PRJ_MIN0 + rnd.PseudoUniform(SGRPROJ_PRJ_MAX0 + 1 -
                                                          SGRPROJ_PRJ_MIN0),
                     SGRPROJ_PRJ_MIN1 + rnd.PseudoUniform(SGRPROJ_PRJ_MAX1 + 1 -
                                                          SGRPROJ_PRJ_MIN1) };
      int eps = rnd.PseudoUniform(1 << SGRPROJ_PARAMS_BITS);
#else
            SVTRandom rnd1(0, SGRPROJ_PRJ_MAX0 + 1 - SGRPROJ_PRJ_MIN0);
            SVTRandom rnd2(0, SGRPROJ_PRJ_MAX1 + 1 - SGRPROJ_PRJ_MIN1);
            SVTRandom eps_rnd(0, 1 << SGRPROJ_PARAMS_BITS);
            int xqd[2] = {SGRPROJ_PRJ_MIN0 + rnd1.random(),
                          SGRPROJ_PRJ_MIN1 + rnd2.random()};
            int eps = eps_rnd.random();
#endif

            // Test various tile sizes around 256x256
            int test_w = max_w - (i / 9);
            int test_h = max_h - (i % 9);

            for (k = 0; k < test_h; k += pu_height)
                for (j = 0; j < test_w; j += pu_width) {
                    int w = AOMMIN(pu_width, test_w - j);
                    int h = AOMMIN(pu_height, test_h - k);
                    uint16_t *input_p = input + k * stride + j;
                    uint16_t *output_p = output + k * out_stride + j;
                    uint16_t *output2_p = output2 + k * out_stride + j;
                    tst_fun_(CONVERT_TO_BYTEPTR(input_p),
                             w,
                             h,
                             stride,
                             eps,
                             xqd,
                             CONVERT_TO_BYTEPTR(output_p),
                             out_stride,
                             tmpbuf,
                             bit_depth,
                             1);
                    apply_selfguided_restoration_c(
                        CONVERT_TO_BYTEPTR(input_p),
                        w,
                        h,
                        stride,
                        eps,
                        xqd,
                        CONVERT_TO_BYTEPTR(output2_p),
                        out_stride,
                        tmpbuf,
                        bit_depth,
                        1);
                }

            for (j = 0; j < test_h; ++j)
                for (k = 0; k < test_w; ++k)
                    ASSERT_EQ(output[j * out_stride + k],
                              output2[j * out_stride + k]);
        }

        aom_free(input_);
        aom_free(output_);
        aom_free(output2_);
        aom_free(tmpbuf);
    }

  private:
    SgrFunc tst_fun_;
};

TEST_P(AV1HighbdSelfguidedFilterTest, CorrectnessTest) {
    RunCorrectnessTest();
}

const int highbd_params_avx2[] = {8, 10, 12};
INSTANTIATE_TEST_CASE_P(
    AVX2, AV1HighbdSelfguidedFilterTest,
    ::testing::Combine(::testing::Values(apply_selfguided_restoration_avx2),
                       ::testing::ValuesIn(highbd_params_avx2)));
}  // namespace
