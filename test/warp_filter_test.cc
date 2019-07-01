/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

/******************************************************************************
 * @file warp_filter_test.cc
 *
 * @brief Unit test of warp filter:
 * TODO: waiting for asm finished
 * - av1_warp_affine_asm
 * - av1_highbd_warp_affine_asm
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbWarpedMotion.h"
#include "EbUnitTestUtility.h"
#include "random.h"
#include "util.h"

namespace {
using std::make_tuple;
using svt_av1_test_tool::SVTRandom;

extern "C" void av1_warp_affine_c(const int32_t* mat, const uint8_t* ref,
                                  int width, int height, int stride,
                                  uint8_t* pred, int p_col, int p_row,
                                  int p_width, int p_height, int p_stride,
                                  int subsampling_x, int subsampling_y,
                                  ConvolveParams* conv_params, int16_t alpha,
                                  int16_t beta, int16_t gamma, int16_t delta);

extern "C" void av1_highbd_warp_affine_c(
    const int32_t* mat, const uint16_t* ref, int width, int height, int stride,
    uint16_t* pred, int p_col, int p_row, int p_width, int p_height,
    int p_stride, int subsampling_x, int subsampling_y, int bd,
    ConvolveParams* conv_params, int16_t alpha, int16_t beta, int16_t gamma,
    int16_t delta);

static const int quant_dist_lookup_table[2][4][2] = {
    {{9, 7}, {11, 5}, {12, 4}, {13, 3}},
    {{7, 9}, {5, 11}, {4, 12}, {3, 13}},
};

class SVTRandomHelper {
  protected:
    SVTRandomHelper() : rnd8_(8, false), rnd16_(16, false) {
    }

    int32_t random_warped_param(int bits) {
        // 1 in 8 chance of generating zero (arbitrarily chosen)
        if (((rnd8_.random()) & 7) == 0)
            return 0;
        // Otherwise, enerate uniform values in the range
        // [-(1 << bits), 1] U [1, 1<<bits]
        int32_t v = 1 + (rnd16_.random() & ((1 << bits) - 1));
        if ((rnd8_.random()) & 1)
            return -v;
        return v;
    }

    void generate_warped_model(int32_t* mat, int16_t* alpha, int16_t* beta,
                               int16_t* gamma, int16_t* delta) {
        while (1) {
            mat[0] = random_warped_param(WARPEDMODEL_PREC_BITS + 6);
            mat[1] = random_warped_param(WARPEDMODEL_PREC_BITS + 6);
            mat[2] = (random_warped_param(WARPEDMODEL_PREC_BITS - 3)) +
                     (1 << WARPEDMODEL_PREC_BITS);
            mat[3] = random_warped_param(WARPEDMODEL_PREC_BITS - 3);
            // 50/50 chance of generating ROTZOOM vs. AFFINE models
            if (rnd8_.random() & 1) {
                // AFFINE
                mat[4] = random_warped_param(WARPEDMODEL_PREC_BITS - 3);
                mat[5] = (random_warped_param(WARPEDMODEL_PREC_BITS - 3)) +
                         (1 << WARPEDMODEL_PREC_BITS);
            } else {
                mat[4] = -mat[3];
                mat[5] = mat[2];
            }

            // Calculate the derived parameters and check that they are suitable
            // for the warp filter.
            assert(mat[2] != 0);

            *alpha = (int16_t)clamp(
                mat[2] - (1 << WARPEDMODEL_PREC_BITS), INT16_MIN, INT16_MAX);
            *beta = (int16_t)clamp(mat[3], INT16_MIN, INT16_MAX);
            *gamma = (int16_t)clamp(
                ((int64_t)mat[4] * (1 << WARPEDMODEL_PREC_BITS)) / mat[2],
                INT16_MIN,
                INT16_MAX);
            *delta = (int16_t)clamp(
                mat[5] - (((int64_t)mat[3] * mat[4] + (mat[2] / 2)) / mat[2]) -
                    (1 << WARPEDMODEL_PREC_BITS),
                INT16_MIN,
                INT16_MAX);

            if ((4 * abs(*alpha) + 7 * abs(*beta) >=
                 (1 << WARPEDMODEL_PREC_BITS)) ||
                (4 * abs(*gamma) + 4 * abs(*delta) >=
                 (1 << WARPEDMODEL_PREC_BITS)))
                continue;

            *alpha = ROUND_POWER_OF_TWO_SIGNED(*alpha, WARP_PARAM_REDUCE_BITS) *
                     (1 << WARP_PARAM_REDUCE_BITS);
            *beta = ROUND_POWER_OF_TWO_SIGNED(*beta, WARP_PARAM_REDUCE_BITS) *
                    (1 << WARP_PARAM_REDUCE_BITS);
            *gamma = ROUND_POWER_OF_TWO_SIGNED(*gamma, WARP_PARAM_REDUCE_BITS) *
                     (1 << WARP_PARAM_REDUCE_BITS);
            *delta = ROUND_POWER_OF_TWO_SIGNED(*delta, WARP_PARAM_REDUCE_BITS) *
                     (1 << WARP_PARAM_REDUCE_BITS);

            // We have a valid model, so finish
            return;
        }
    }

  protected:
    SVTRandom rnd8_;
    SVTRandom rnd16_;
};

typedef void (*warp_affine_func)(const int32_t* mat, const uint8_t* ref,
                                 int width, int height, int stride,
                                 uint8_t* pred, int p_col, int p_row,
                                 int p_width, int p_height, int p_stride,
                                 int subsampling_x, int subsampling_y,
                                 ConvolveParams* conv_params, int16_t alpha,
                                 int16_t beta, int16_t gamma, int16_t delta);

typedef std::tuple<int, int, int, warp_affine_func, int /* constant 8 */>
    WarpTestParam;

typedef void (*highbd_warp_affine_func)(const int32_t* mat, const uint16_t* ref,
                                        int width, int height, int stride,
                                        uint16_t* pred, int p_col, int p_row,
                                        int p_width, int p_height, int p_stride,
                                        int subsampling_x, int subsampling_y,
                                        int bd, ConvolveParams* conv_params,
                                        int16_t alpha, int16_t beta,
                                        int16_t gamma, int16_t delta);

typedef std::tuple<int, int, int, highbd_warp_affine_func, int>
    HighbdWarpTestParam;

static const int w = 128, h = 128;
static const int border = 16;
static const int stride = w + 2 * border;

template <typename Sample, typename FuncType, typename ParamType>
class AV1WarpFilterTestTemplate : public ::testing::TestWithParam<ParamType>,
                                  public SVTRandomHelper {
  protected:
    void SetUp() override {
        func_tst_ = TEST_GET_PARAM(3);
        bd_ = TEST_GET_PARAM(4);
        input_ = nullptr;
        output_tst_ = nullptr;
        output_ref_ = nullptr;
        out_width_ = TEST_GET_PARAM(0);
        out_height_ = TEST_GET_PARAM(1);
        test_times_ = TEST_GET_PARAM(2);
        rnd8_.reset();
        rnd16_.reset();

        // The warp functions always write rows with widths that are multiples
        // of 8. So to avoid a buffer overflow, we may need to pad rows to a
        // multiple of 8.
        int output_n = ((out_width_ + 7) & ~7) * out_height_;
        input_ = new Sample[h * stride];
        output_tst_ = new Sample[output_n];
        output_ref_ = new Sample[output_n];
        for (int i = 0; i < output_n; ++i)
            output_tst_[i] = output_ref_[i] = gen_random();
    }

    void TearDown() override {
        if (input_) {
            delete[] input_;
            input_ = nullptr;
        }
        if (output_tst_) {
            delete[] output_tst_;
            output_tst_ = nullptr;
        }
        if (output_ref_) {
            delete[] output_ref_;
            output_ref_ = nullptr;
        }
        aom_clear_system_state();
    }

    template <typename T>
    void compare_outputs(const T* output1, const T* output2, const int width,
                         const int height, const int it) {
        for (int i = 0; i < width * height; ++i)
            ASSERT_EQ(output1[i], output2[i])
                << "Pixel mismatch at index " << i << " = (" << (i % width)
                << ", " << (i / height) << ") on iteration " << it;
    }

    void prepare_data() {
        const int mask = (1 << bd_) - 1;
        Sample* input = input_ + border;
        // Generate an input block and extend its borders horizontally
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
                input[r * stride + c] = gen_random() & mask;
        for (int r = 0; r < h; ++r) {
            for (int c = 0; c < border; ++c) {
                input[r * stride - border + c] = input[r * stride];
                input[r * stride + w + c] = input[r * stride + (w - 1)];
            }
        }
    }

    virtual Sample gen_random() = 0;
    virtual void call_warp_affine_test(const int sub_x, const int sub_y,
                                       ConvolveParams& conv_params,
                                       const int32_t mat[], const int16_t alpha,
                                       const int16_t beta, const int16_t gamma,
                                       const int16_t delta) = 0;
    virtual void call_warp_affine_reference(
        const int sub_x, const int sub_y, ConvolveParams& conv_params,
        const int32_t mat[], const int16_t alpha, const int16_t beta,
        const int16_t gamma, const int16_t delta) = 0;

    void do_run_and_compare(bool use_no_round, int sub_x, int sub_y,
                            uint16_t* dsta, uint16_t* dstb, const int times) {
        int32_t mat[8];
        int16_t alpha, beta, gamma, delta;
        generate_warped_model(mat, &alpha, &beta, &gamma, &delta);

        ConvolveParams conv_params;
        for (int ii = 0; ii < 2; ++ii) {
            for (int jj = 0; jj < 5; ++jj) {
                for (int do_average = 0; do_average <= 1; ++do_average) {
                    if (use_no_round) {
                        conv_params = get_conv_params_no_round(
                            0, do_average, 0, dsta, out_width_, 1, bd_);
                    } else {
                        conv_params = get_conv_params(0, 0, 0, bd_);
                    }
                    if (jj >= 4) {
                        conv_params.use_jnt_comp_avg = 0;
                    } else {
                        conv_params.use_jnt_comp_avg = 1;
                        conv_params.fwd_offset =
                            quant_dist_lookup_table[ii][jj][0];
                        conv_params.bck_offset =
                            quant_dist_lookup_table[ii][jj][1];
                    }
                    call_warp_affine_test(sub_x,
                                          sub_y,
                                          conv_params,
                                          mat,
                                          alpha,
                                          beta,
                                          gamma,
                                          delta);

                    if (use_no_round) {
                        conv_params = get_conv_params_no_round(
                            0, do_average, 0, dstb, out_width_, 1, bd_);
                    }
                    if (jj >= 4) {
                        conv_params.use_jnt_comp_avg = 0;
                    } else {
                        conv_params.use_jnt_comp_avg = 1;
                        conv_params.fwd_offset =
                            quant_dist_lookup_table[ii][jj][0];
                        conv_params.bck_offset =
                            quant_dist_lookup_table[ii][jj][1];
                    }
                    call_warp_affine_reference(sub_x,
                                               sub_y,
                                               conv_params,
                                               mat,
                                               alpha,
                                               beta,
                                               gamma,
                                               delta);

                    if (use_no_round) {
                        compare_outputs(
                            dsta, dstb, out_width_, out_height_, times);
                        compare_outputs(output_tst_,
                                        output_ref_,
                                        out_width_,
                                        out_height_,
                                        times);
                    } else {
                        compare_outputs(output_tst_,
                                        output_ref_,
                                        out_width_,
                                        out_height_,
                                        times);
                    }
                }
            }
        }
    }

    virtual void run_check_output() {
        int output_n = ((out_width_ + 7) & ~7) * out_height_;
        uint16_t* dsta = new uint16_t[output_n];
        uint16_t* dstb = new uint16_t[output_n];
        ASSERT_NE(dsta, nullptr);
        ASSERT_NE(dstb, nullptr);
        for (int i = 0; i < test_times_; ++i) {
            prepare_data();
            const int use_no_round = rnd8_.random() & 1;
            for (int sub_x = 0; sub_x < 2; ++sub_x) {
                for (int sub_y = 0; sub_y < 2; ++sub_y)
                    do_run_and_compare(
                        use_no_round, sub_x, sub_y, dsta, dstb, i);
            }
        }

        delete[] dsta;
        delete[] dstb;
    }

  protected:
    FuncType func_tst_;  /**< function to test */
    int bd_;             /**< bit-depth */
    Sample* input_;      /**< buffer of input data */
    Sample* output_tst_; /**< buffer of test output data */
    Sample* output_ref_; /**< buffer of reference output data */
    int out_width_;      /**< width of output */
    int out_height_;     /**< heihgt of output */
    int test_times_;     /**< total times to test */
};

class AV1WarpFilterTest
    : public AV1WarpFilterTestTemplate<uint8_t, warp_affine_func,
                                       WarpTestParam> {
  public:
    static ::testing::internal::ParamGenerator<WarpTestParam> BuildParams(
        warp_affine_func filter) {
        const WarpTestParam params[] = {
            make_tuple(4, 4, 500, filter, 8),
            make_tuple(8, 8, 500, filter, 8),
            make_tuple(64, 64, 100, filter, 8),
            make_tuple(4, 16, 200, filter, 8),
            make_tuple(32, 8, 100, filter, 8),
        };
        return ::testing::ValuesIn(params);
    }

  protected:
    uint8_t gen_random() override {
        return (uint8_t)rnd8_.random();
    }
    void call_warp_affine_test(const int sub_x, const int sub_y,
                               ConvolveParams& conv_params, const int32_t mat[],
                               const int16_t alpha, const int16_t beta,
                               const int16_t gamma,
                               const int16_t delta) override {
        uint8_t* input = input_ + border;
        func_tst_(mat,
                  input,
                  w,
                  h,
                  stride,
                  output_tst_,
                  32,
                  32,
                  out_width_,
                  out_height_,
                  out_width_,
                  sub_x,
                  sub_y,
                  &conv_params,
                  alpha,
                  beta,
                  gamma,
                  delta);
    }
    void call_warp_affine_reference(const int sub_x, const int sub_y,
                                    ConvolveParams& conv_params,
                                    const int32_t mat[], const int16_t alpha,
                                    const int16_t beta, const int16_t gamma,
                                    const int16_t delta) override {
        uint8_t* input = input_ + border;
        av1_warp_affine_c(mat,
                          input,
                          w,
                          h,
                          stride,
                          output_ref_,
                          32,
                          32,
                          out_width_,
                          out_height_,
                          out_width_,
                          sub_x,
                          sub_y,
                          &conv_params,
                          alpha,
                          beta,
                          gamma,
                          delta);
    }
};

class AV1HBDWarpFilterTest
    : public AV1WarpFilterTestTemplate<uint16_t, highbd_warp_affine_func,
                                       HighbdWarpTestParam> {
  public:
    static ::testing::internal::ParamGenerator<HighbdWarpTestParam> BuildParams(
        highbd_warp_affine_func filter) {
        const HighbdWarpTestParam params[] = {
            make_tuple(4, 4, 100, filter, 8),
            make_tuple(8, 8, 100, filter, 8),
            make_tuple(64, 64, 100, filter, 8),
            make_tuple(4, 16, 100, filter, 8),
            make_tuple(32, 8, 100, filter, 8),
            make_tuple(4, 4, 100, filter, 10),
            make_tuple(8, 8, 100, filter, 10),
            make_tuple(64, 64, 100, filter, 10),
            make_tuple(4, 16, 100, filter, 10),
            make_tuple(32, 8, 100, filter, 10),
            make_tuple(4, 4, 100, filter, 12),
            make_tuple(8, 8, 100, filter, 12),
            make_tuple(64, 64, 100, filter, 12),
            make_tuple(4, 16, 100, filter, 12),
            make_tuple(32, 8, 100, filter, 12),
        };
        return ::testing::ValuesIn(params);
    }

  protected:
    uint16_t gen_random() override {
        return (uint16_t)rnd16_.random();
    }
    void call_warp_affine_test(const int sub_x, const int sub_y,
                               ConvolveParams& conv_params, const int32_t mat[],
                               const int16_t alpha, const int16_t beta,
                               const int16_t gamma,
                               const int16_t delta) override {
        uint16_t* input = input_ + border;
        func_tst_(mat,
                  input,
                  w,
                  h,
                  stride,
                  output_tst_,
                  32,
                  32,
                  out_width_,
                  out_height_,
                  out_width_,
                  sub_x,
                  sub_y,
                  bd_,
                  &conv_params,
                  alpha,
                  beta,
                  gamma,
                  delta);
    }
    void call_warp_affine_reference(const int sub_x, const int sub_y,
                                    ConvolveParams& conv_params,
                                    const int32_t mat[], const int16_t alpha,
                                    const int16_t beta, const int16_t gamma,
                                    const int16_t delta) override {
        uint16_t* input = input_ + border;
        av1_highbd_warp_affine_c(mat,
                                 input,
                                 w,
                                 h,
                                 stride,
                                 output_ref_,
                                 32,
                                 32,
                                 out_width_,
                                 out_height_,
                                 out_width_,
                                 sub_x,
                                 sub_y,
                                 bd_,
                                 &conv_params,
                                 alpha,
                                 beta,
                                 gamma,
                                 delta);
    }
};

TEST_P(AV1WarpFilterTest, run_check_output) {
    run_check_output();
}

INSTANTIATE_TEST_CASE_P(AV1, AV1WarpFilterTest,
                        AV1WarpFilterTest::BuildParams(av1_warp_affine_c));

TEST_P(AV1HBDWarpFilterTest, run_check_output) {
    run_check_output();
}

INSTANTIATE_TEST_CASE_P(
    AV1, AV1HBDWarpFilterTest,
    AV1HBDWarpFilterTest::BuildParams(av1_highbd_warp_affine_c));

}  // namespace
