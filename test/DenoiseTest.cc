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
 * @file DenoiseTest.cc
 *
 * @brief Unit test of wiener demoise 2d function:
 * - aom_wiener_denoise_2d
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "gtest/gtest.h"
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "random.h"
#include "util.h"
#include "noise_util.h"
#include "noise_model.h"
#include <vector>

namespace {
using std::make_tuple;
using svt_av1_test_tool::SVTRandom;

static const int kBlockSize = 32;
static const int kBlockSizeChroma = 16;
static const int kWidth = 256;
static const int kHeight = 256;

/** setup_test_env is implemented in test/TestEnv.c */
extern "C" void setup_test_env(EbAsm asm_type);

/** fclamp is copied from nois_model.c */
static inline double fclamp(double value, double low, double high) {
    return value < low ? low : (value > high ? high : value);
}

/** same data struct as aom_noise_tx_t, renamed to avoid name conflict*/
struct test_noise_tx_t {
    DECLARE_ALIGNED(32, float, *tx_block);
    DECLARE_ALIGNED(32, float, *temp);
    int32_t block_size;
    void (*fft)(const float*, float*, float*);
    void (*ifft)(const float*, float*, float*);
};

/**for aom_noise_tx_add_energy is deleted in noise_util.c, test copied from
 * aom*/
static inline void aom_noise_tx_add_energy(
    const struct test_noise_tx_t* noise_tx, float* psd) {
    const int block_size = noise_tx->block_size;
    for (int yb = 0; yb < block_size; ++yb) {
        for (int xb = 0; xb <= block_size / 2; ++xb) {
            float* c = noise_tx->tx_block + 2 * (yb * block_size + xb);
            psd[yb * block_size + xb] += c[0] * c[0] + c[1] * c[1];
        }
    }
}

typedef std::tuple<int,  /**< bit-depth */
                   EbAsm /**< asm code accelerate of fft/ifft for test */
                   >
    WienerDenoisearam;

template <typename Sample>
class WienerDenoiseTest : public ::testing::TestWithParam<WienerDenoisearam> {
  protected:
    WienerDenoiseTest() : rnd31_(31, false) {
        bd_ = TEST_GET_PARAM(0);
        scale_noise_ = 1 << (bd_ - 8);
        setup_test_env(TEST_GET_PARAM(1));
    }

    // Return normally distrbuted values with standard deviation of sigma.
    double randn(double sigma) {
        while (1) {
            const double u = 2.0 * ((double)rnd31_.random() / (1u << 31)) - 1.0;
            const double v = 2.0 * ((double)rnd31_.random() / (1u << 31)) - 1.0;
            const double s = u * u + v * v;
            if (s > 0 && s < 1) {
                return sigma * (u * sqrt(-2.0 * log(s) / s));
            }
        }
        return 0;
    }

    // Synthesizes noise using the auto-regressive filter of the given lag,
    // with the provided n coefficients sampled at the given coords.
    void noise_synth(int lag, int n, const int (*coords)[2],
                     const double* coeffs, double* data, int w, int h) {
        const int pad_size = 3 * lag;
        const int padded_w = w + pad_size;
        const int padded_h = h + pad_size;
        int x = 0, y = 0;
        std::vector<double> padded(padded_w * padded_h);

        for (y = 0; y < padded_h; ++y) {
            for (x = 0; x < padded_w; ++x) {
                padded[y * padded_w + x] = randn(1.0);
            }
        }
        for (y = lag; y < padded_h; ++y) {
            for (x = lag; x < padded_w; ++x) {
                double sum = 0;
                int i = 0;
                for (i = 0; i < n; ++i) {
                    const int dx = coords[i][0];
                    const int dy = coords[i][1];
                    sum += padded[(y + dy) * padded_w + (x + dx)] * coeffs[i];
                }
                padded[y * padded_w + x] += sum;
            }
        }
        // Copy over the padded rows to the output
        for (y = 0; y < h; ++y) {
            memcpy(data + y * w, &padded[0] + y * padded_w, sizeof(*data) * w);
        }
    }

    std::vector<float> get_noise_psd(double* noise, int width, int height,
                                     int block_size) {
        float* block =
            (float*)aom_memalign(32, block_size * block_size * sizeof(block));
        std::vector<float> psd(block_size * block_size);
        int num_blocks = 0;
        struct aom_noise_tx_t* tx = aom_noise_tx_malloc(block_size);
        for (int y = 0; y <= height - block_size; y += block_size / 2) {
            for (int x = 0; x <= width - block_size; x += block_size / 2) {
                for (int yy = 0; yy < block_size; ++yy) {
                    for (int xx = 0; xx < block_size; ++xx) {
                        block[yy * block_size + xx] =
                            (float)noise[(y + yy) * width + x + xx];
                    }
                }
                aom_noise_tx_forward(tx, &block[0]);
                aom_noise_tx_add_energy(
                    reinterpret_cast<const test_noise_tx_t*>(tx), &psd[0]);
                num_blocks++;
            }
        }
        for (int yy = 0; yy < block_size; ++yy) {
            for (int xx = 0; xx <= block_size / 2; ++xx) {
                psd[yy * block_size + xx] /= num_blocks;
            }
        }
        // Fill in the data that is missing due to symmetries
        for (int xx = 1; xx < block_size / 2; ++xx) {
            psd[(block_size - xx)] = psd[xx];
        }
        for (int yy = 1; yy < block_size; ++yy) {
            for (int xx = 1; xx < block_size / 2; ++xx) {
                psd[(block_size - yy) * block_size + (block_size - xx)] =
                    psd[yy * block_size + xx];
            }
        }
        aom_noise_tx_free(tx);
        aom_free(block);
        return psd;
    }

    void SetUp() override {
        const float kNoiseLevel = 5.f;
        const float kStd = 4.0;
        const double kMaxValue = (1 << bd_) - 1;

        chroma_sub_[0] = 1;
        chroma_sub_[1] = 1;
        stride_[0] = kWidth;
        stride_[1] = kWidth / 2;
        stride_[2] = kWidth / 2;
        for (int k = 0; k < 3; ++k) {
            data_[k].resize(kWidth * kHeight);
            denoised_[k].resize(kWidth * kHeight);
            noise_psd_[k].resize(kBlockSize * kBlockSize);
        }

        const double kCoeffsY[] = {0.0406,
                                   -0.116,
                                   -0.078,
                                   -0.152,
                                   0.0033,
                                   -0.093,
                                   0.048,
                                   0.404,
                                   0.2353,
                                   -0.035,
                                   -0.093,
                                   0.441};
        const int kCoords[12][2] = {{-2, -2},
                                    {-1, -2},
                                    {0, -2},
                                    {1, -2},
                                    {2, -2},
                                    {-2, -1},
                                    {-1, -1},
                                    {0, -1},
                                    {1, -1},
                                    {2, -1},
                                    {-2, 0},
                                    {-1, 0}};
        const int kLag = 2;
        const int kLength = 12;

        std::vector<double> noise(kWidth * kHeight);
        noise_synth(
            kLag, kLength, kCoords, kCoeffsY, &noise[0], kWidth, kHeight);
        noise_psd_[0] = get_noise_psd(&noise[0], kWidth, kHeight, kBlockSize);
        for (int i = 0; i < kBlockSize * kBlockSize; ++i) {
            noise_psd_[0][i] =
                (float)(noise_psd_[0][i] * kStd * kStd * scale_noise_ *
                        scale_noise_ / (kMaxValue * kMaxValue));
        }

        float psd_value =
            aom_noise_psd_get_default_value(kBlockSizeChroma, kNoiseLevel);
        for (int i = 0; i < kBlockSizeChroma * kBlockSizeChroma; ++i) {
            noise_psd_[1][i] = psd_value;
            noise_psd_[2][i] = psd_value;
        }
        for (int y = 0; y < kHeight; ++y) {
            for (int x = 0; x < kWidth; ++x) {
                data_[0][y * stride_[0] + x] = (Sample)fclamp(
                    (x + noise[y * stride_[0] + x] * kStd) * scale_noise_,
                    0,
                    kMaxValue);
            }
        }

        for (int c = 1; c < 3; ++c) {
            for (int y = 0; y < (kHeight >> 1); ++y) {
                for (int x = 0; x < (kWidth >> 1); ++x) {
                    data_[c][y * stride_[c] + x] = (Sample)fclamp(
                        (x + randn(kStd)) * scale_noise_, 0, kMaxValue);
                }
            }
        }
        for (int k = 0; k < 3; ++k) {
            noise_psd_ptrs_[k] = &noise_psd_[k][0];
        }
    }

  protected:
    void run_invalid_block_size() {
        const uint8_t* const data_ptrs[3] = {
            reinterpret_cast<uint8_t*>(&data_[0][0]),
            reinterpret_cast<uint8_t*>(&data_[1][0]),
            reinterpret_cast<uint8_t*>(&data_[2][0]),
        };
        uint8_t* denoised_ptrs[3] = {
            reinterpret_cast<uint8_t*>(&denoised_[0][0]),
            reinterpret_cast<uint8_t*>(&denoised_[1][0]),
            reinterpret_cast<uint8_t*>(&denoised_[2][0]),
        };
        EXPECT_EQ(aom_wiener_denoise_2d(data_ptrs,
                                        denoised_ptrs,
                                        kWidth,
                                        kHeight,
                                        stride_,
                                        chroma_sub_,
                                        noise_psd_ptrs_,
                                        18,
                                        bd_,
                                        use_hdb_),
                  0);
        EXPECT_EQ(aom_wiener_denoise_2d(data_ptrs,
                                        denoised_ptrs,
                                        kWidth,
                                        kHeight,
                                        stride_,
                                        chroma_sub_,
                                        noise_psd_ptrs_,
                                        48,
                                        bd_,
                                        use_hdb_),
                  0);
        EXPECT_EQ(aom_wiener_denoise_2d(data_ptrs,
                                        denoised_ptrs,
                                        kWidth,
                                        kHeight,
                                        stride_,
                                        chroma_sub_,
                                        noise_psd_ptrs_,
                                        64,
                                        bd_,
                                        use_hdb_),
                  0);
    }

    void run_invalid_chroma_subsampling() {
        const uint8_t* const data_ptrs[3] = {
            reinterpret_cast<uint8_t*>(&data_[0][0]),
            reinterpret_cast<uint8_t*>(&data_[1][0]),
            reinterpret_cast<uint8_t*>(&data_[2][0]),
        };
        uint8_t* denoised_ptrs[3] = {
            reinterpret_cast<uint8_t*>(&denoised_[0][0]),
            reinterpret_cast<uint8_t*>(&denoised_[1][0]),
            reinterpret_cast<uint8_t*>(&denoised_[2][0]),
        };
        int chroma_sub[2] = {1, 0};
        EXPECT_EQ(aom_wiener_denoise_2d(data_ptrs,
                                        denoised_ptrs,
                                        kWidth,
                                        kHeight,
                                        stride_,
                                        chroma_sub,
                                        noise_psd_ptrs_,
                                        32,
                                        bd_,
                                        use_hdb_),
                  0);

        chroma_sub[0] = 0;
        chroma_sub[1] = 1;
        EXPECT_EQ(aom_wiener_denoise_2d(data_ptrs,
                                        denoised_ptrs,
                                        kWidth,
                                        kHeight,
                                        stride_,
                                        chroma_sub,
                                        noise_psd_ptrs_,
                                        32,
                                        bd_,
                                        use_hdb_),
                  0);
    }

    void run_gradient_test() {
        const uint8_t* const data_ptrs[3] = {
            reinterpret_cast<uint8_t*>(&data_[0][0]),
            reinterpret_cast<uint8_t*>(&data_[1][0]),
            reinterpret_cast<uint8_t*>(&data_[2][0]),
        };
        uint8_t* denoised_ptrs[3] = {
            reinterpret_cast<uint8_t*>(&denoised_[0][0]),
            reinterpret_cast<uint8_t*>(&denoised_[1][0]),
            reinterpret_cast<uint8_t*>(&denoised_[2][0]),
        };
        const int ret = aom_wiener_denoise_2d(data_ptrs,
                                              denoised_ptrs,
                                              kWidth,
                                              kHeight,
                                              stride_,
                                              chroma_sub_,
                                              noise_psd_ptrs_,
                                              kBlockSize,
                                              bd_,
                                              use_hdb_);
        EXPECT_EQ(ret, 1);

        // Check the noise on the denoised image (from the analytical gradient)
        // and make sure that it is less than what we added.
        for (int c = 0; c < 3; ++c) {
            std::vector<double> measured_noise(kWidth * kHeight);

            double var = 0;
            const int shift = (c > 0);
            for (int x = 0; x < (kWidth >> shift); ++x) {
                for (int y = 0; y < (kHeight >> shift); ++y) {
                    const double diff =
                        denoised_[c][y * stride_[c] + x] - x * scale_noise_;
                    var += diff * diff;
                    measured_noise[y * kWidth + x] = diff;
                }
            }
            var /= (kWidth * kHeight);
            const double std = sqrt(max(0.0, var));
            EXPECT_LE(std, 1.25f * scale_noise_);
            if (c == 0) {
                std::vector<float> measured_psd = get_noise_psd(
                    &measured_noise[0], kWidth, kHeight, kBlockSize);
                std::vector<double> measured_psd_d(kBlockSize * kBlockSize);
                std::vector<double> noise_psd_d(kBlockSize * kBlockSize);
                std::copy(measured_psd.begin(),
                          measured_psd.end(),
                          measured_psd_d.begin());
                std::copy(noise_psd_[0].begin(),
                          noise_psd_[0].end(),
                          noise_psd_d.begin());
                EXPECT_LT(
                    aom_normalized_cross_correlation(&measured_psd_d[0],
                                                     &noise_psd_d[0],
                                                     (int)(noise_psd_d.size())),
                    0.35);
            }
        }
    }

  protected:
    std::vector<Sample> data_[3];
    std::vector<Sample> denoised_[3];
    std::vector<float> noise_psd_[3];
    int chroma_sub_[2];
    float* noise_psd_ptrs_[3];
    int stride_[3];

    SVTRandom rnd31_; /**< random tool */
    uint8_t bd_;      /**< bit depth */
    int32_t use_hdb_; /**< flag of use HBD or not */
    int scale_noise_; /**< scale factor of noise */
};

class WienerDenoiseLbdTest : public WienerDenoiseTest<uint8_t> {
  protected:
    WienerDenoiseLbdTest() {
        use_hdb_ = 0;
    }
};

TEST_P(WienerDenoiseLbdTest, run_invalid_block_size) {
    run_invalid_block_size();
}
TEST_P(WienerDenoiseLbdTest, run_invalid_chroma_subsampling) {
    run_invalid_chroma_subsampling();
}
TEST_P(WienerDenoiseLbdTest, run_gradient_test) {
    run_gradient_test();
}

static const WienerDenoisearam lbd_params[] = {make_tuple(8, ASM_NON_AVX2),
                                               make_tuple(8, ASM_AVX2),
                                               make_tuple(8, ASM_TYPE_TOTAL)};

INSTANTIATE_TEST_CASE_P(Denoise, WienerDenoiseLbdTest,
                        ::testing::ValuesIn(lbd_params));

class WienerDenoiseHbdTest : public WienerDenoiseTest<uint16_t> {
  protected:
    WienerDenoiseHbdTest() {
        use_hdb_ = 1;
    }
};

TEST_P(WienerDenoiseHbdTest, run_invalid_block_size) {
    run_invalid_block_size();
}
TEST_P(WienerDenoiseHbdTest, run_invalid_chroma_subsampling) {
    run_invalid_chroma_subsampling();
}
TEST_P(WienerDenoiseHbdTest, run_gradient_test) {
    run_gradient_test();
}

static const WienerDenoisearam hbd_params[] = {make_tuple(8, ASM_NON_AVX2),
                                               make_tuple(8, ASM_AVX2),
                                               make_tuple(8, ASM_TYPE_TOTAL),
                                               make_tuple(10, ASM_NON_AVX2),
                                               make_tuple(10, ASM_AVX2),
                                               make_tuple(10, ASM_TYPE_TOTAL),
                                               make_tuple(12, ASM_NON_AVX2),
                                               make_tuple(12, ASM_AVX2),
                                               make_tuple(12, ASM_TYPE_TOTAL)};

INSTANTIATE_TEST_CASE_P(Denoise, WienerDenoiseHbdTest,
                        ::testing::ValuesIn(hbd_params));

}  // namespace
