/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

/******************************************************************************
 * @file ObmcVarianceTest.cc
 *
 * @brief Unit test for OBMC sad
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

namespace {
#define MAX_BLOCK_SIZE (128 * 256)
#define MAX_SB_SQUARE (MAX_BLOCK_SIZE * MAX_BLOCK_SIZE)

static const int kMaskMax = 64;
static const int kIterations = 1000;

using svt_av1_test_tool::SVTRandom;

typedef uint32_t (*EB_OBMC_VARIANCE_FCN)(const uint8_t* pre, int pre_stride,
                                         const int32_t* wsrc,
                                         const int32_t* mask, uint32_t* sse);

typedef std::tuple<EB_OBMC_VARIANCE_FCN, EB_OBMC_VARIANCE_FCN>
    ObmcVarianceTestParam;

class ObmcVarianceTest
    : public ::testing::TestWithParam<ObmcVarianceTestParam> {
  protected:
    void check_obmc_random() {
        EB_OBMC_VARIANCE_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_VARIANCE_FCN tst_fcn = TEST_GET_PARAM(1);
        ASSERT_NE(ref_fcn, nullptr) << "can not find reference function";
        ASSERT_NE(tst_fcn, nullptr) << "can not find test function";

        uint8_t* pre =
            (uint8_t*)aom_memalign(32, sizeof(uint8_t) * MAX_SB_SQUARE);
        int32_t* wsrc =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);
        int32_t* mask =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);

        SVTRandom rnd(8, false);
        SVTRandom rnd_stride(0, MAX_SB_SIZE + 1);
        SVTRandom rnd_mask(0, kMaskMax * kMaskMax + 1);
        for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
            for (int i = 0; i < MAX_SB_SQUARE; ++i) {
                pre[i] = rnd.random();
                wsrc[i] = rnd.random() * rnd_mask.random();
                mask[i] = rnd_mask.random();
            }

            const int pre_stride = rnd_stride.random();
            uint32_t ref_sse, tst_sse;
            const uint32_t ref_res =
                ref_fcn(pre, pre_stride, wsrc, mask, &ref_sse);
            const uint32_t tst_res =
                tst_fcn(pre, pre_stride, wsrc, mask, &tst_sse);
            ASSERT_EQ(tst_res, ref_res);
            ASSERT_EQ(tst_sse, ref_sse);
        }
        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }

    void check_obmc_exterme_values() {
        EB_OBMC_VARIANCE_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_VARIANCE_FCN tst_fcn = TEST_GET_PARAM(1);
        ASSERT_NE(ref_fcn, nullptr) << "can not find reference function";
        ASSERT_NE(tst_fcn, nullptr) << "can not find test function";

        uint8_t* pre =
            (uint8_t*)aom_memalign(32, sizeof(uint8_t) * MAX_SB_SQUARE);
        int32_t* wsrc =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);
        int32_t* mask =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);

        for (int iter = 0; iter < MAX_SB_SIZE && !HasFatalFailure(); ++iter) {
            for (int i = 0; i < MAX_SB_SQUARE; ++i) {
                pre[i] = UINT8_MAX;
                wsrc[i] = UINT8_MAX * kMaskMax * kMaskMax;
                mask[i] = kMaskMax * kMaskMax;
            }

            const int pre_stride = iter;
            uint32_t ref_sse, tst_sse;
            const uint32_t ref_res =
                ref_fcn(pre, pre_stride, wsrc, mask, &ref_sse);
            const uint32_t tst_res =
                tst_fcn(pre, pre_stride, wsrc, mask, &tst_sse);
            ASSERT_EQ(tst_res, ref_res);
            ASSERT_EQ(tst_sse, ref_sse);
        }
        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }
};

TEST_P(ObmcVarianceTest, RandomValues) {
    check_obmc_random();
}

TEST_P(ObmcVarianceTest, ExtremeValues) {
    check_obmc_exterme_values();
}

/* TODO: fill the target function and reference function in test vector
const ObmcVarianceTest::ParamType sse4_functions[] = {
  TestFuncs(aom_obmc_variance128x128_c, aom_obmc_variance128x128_sse4_1),
  TestFuncs(aom_obmc_variance128x64_c, aom_obmc_variance128x64_sse4_1),
  TestFuncs(aom_obmc_variance64x128_c, aom_obmc_variance64x128_sse4_1),
  TestFuncs(aom_obmc_variance64x64_c, aom_obmc_variance64x64_sse4_1),
  TestFuncs(aom_obmc_variance64x32_c, aom_obmc_variance64x32_sse4_1),
  TestFuncs(aom_obmc_variance32x64_c, aom_obmc_variance32x64_sse4_1),
  TestFuncs(aom_obmc_variance32x32_c, aom_obmc_variance32x32_sse4_1),
  TestFuncs(aom_obmc_variance32x16_c, aom_obmc_variance32x16_sse4_1),
  TestFuncs(aom_obmc_variance16x32_c, aom_obmc_variance16x32_sse4_1),
  TestFuncs(aom_obmc_variance16x16_c, aom_obmc_variance16x16_sse4_1),
  TestFuncs(aom_obmc_variance16x8_c, aom_obmc_variance16x8_sse4_1),
  TestFuncs(aom_obmc_variance8x16_c, aom_obmc_variance8x16_sse4_1),
  TestFuncs(aom_obmc_variance8x8_c, aom_obmc_variance8x8_sse4_1),
  TestFuncs(aom_obmc_variance8x4_c, aom_obmc_variance8x4_sse4_1),
  TestFuncs(aom_obmc_variance4x8_c, aom_obmc_variance4x8_sse4_1),
  TestFuncs(aom_obmc_variance4x4_c, aom_obmc_variance4x4_sse4_1)
};
*/
ObmcVarianceTestParam TEST_PARAM[] = {
    ObmcVarianceTestParam(nullptr, nullptr),
};

INSTANTIATE_TEST_CASE_P(OBMC, ObmcVarianceTest,
                        ::testing::ValuesIn(TEST_PARAM));

////////////////////////////////////////////////////////////////////////////////
// High bit-depth
////////////////////////////////////////////////////////////////////////////////

typedef std::tuple<EB_OBMC_VARIANCE_FCN, EB_OBMC_VARIANCE_FCN, int>
    HBDObmcVarianceTestParam;

class HBDObmcVarianceTest
    : public ::testing::TestWithParam<HBDObmcVarianceTestParam> {
  protected:
    void check_highbd_obmc_random() {
        EB_OBMC_VARIANCE_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_VARIANCE_FCN tst_fcn = TEST_GET_PARAM(1);
        ASSERT_NE(ref_fcn, nullptr) << "can not find reference function";
        ASSERT_NE(tst_fcn, nullptr) << "can not find test function";

        uint16_t* pre =
            (uint16_t*)aom_memalign(32, sizeof(uint16_t) * MAX_SB_SQUARE);
        int32_t* wsrc =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);
        int32_t* mask =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);

        int bd = TEST_GET_PARAM(2);
        SVTRandom rnd(bd, false);
        SVTRandom rnd_stride(0, MAX_SB_SIZE + 1);
        SVTRandom rnd_mask(0, kMaskMax * kMaskMax + 1);
        for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
            for (int i = 0; i < MAX_SB_SQUARE; ++i) {
                pre[i] = rnd.random();
                wsrc[i] = rnd.random() * rnd_mask.random();
                mask[i] = rnd_mask.random();
            }

            const int pre_stride = rnd_stride.random();
            uint32_t ref_sse, tst_sse;
            const uint32_t ref_res =
                ref_fcn((uint8_t*)pre, pre_stride, wsrc, mask, &ref_sse);
            const uint32_t tst_res =
                tst_fcn((uint8_t*)pre, pre_stride, wsrc, mask, &tst_sse);
            ASSERT_EQ(tst_res, ref_res);
            ASSERT_EQ(tst_sse, ref_sse);
        }
        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }

    void check_highbd_obmc_exterme_values() {
        EB_OBMC_VARIANCE_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_VARIANCE_FCN tst_fcn = TEST_GET_PARAM(1);
        ASSERT_NE(ref_fcn, nullptr) << "can not find reference function";
        ASSERT_NE(tst_fcn, nullptr) << "can not find test function";

        uint16_t* pre =
            (uint16_t*)aom_memalign(32, sizeof(uint16_t) * MAX_SB_SQUARE);
        int32_t* wsrc =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);
        int32_t* mask =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);

        int bd = TEST_GET_PARAM(2);
        for (int iter = 0; iter < MAX_SB_SIZE && !HasFatalFailure(); ++iter) {
            for (int i = 0; i < MAX_SB_SQUARE; ++i) {
                pre[i] = (1 << bd) - 1;
                wsrc[i] = ((1 << bd) - 1) * kMaskMax * kMaskMax;
                mask[i] = kMaskMax * kMaskMax;
            }

            const int pre_stride = iter;
            uint32_t ref_sse, tst_sse;
            const uint32_t ref_res =
                ref_fcn((uint8_t*)pre, pre_stride, wsrc, mask, &ref_sse);
            const uint32_t tst_res =
                tst_fcn((uint8_t*)pre, pre_stride, wsrc, mask, &tst_sse);
            ASSERT_EQ(tst_res, ref_res);
            ASSERT_EQ(tst_sse, ref_sse);
        }

        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }
};

TEST_P(HBDObmcVarianceTest, RandomValues) {
    check_highbd_obmc_random();
}

TEST_P(HBDObmcVarianceTest, ExtremeValues) {
    check_highbd_obmc_exterme_values();
}

/* TODO: fill the target function and reference function in test vector
ObmcVarianceHBDTest::ParamType sse4_functions_hbd[] = {
  TestFuncs(aom_highbd_obmc_variance128x128_c,
                        aom_highbd_obmc_variance128x128_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance128x64_c,
                        aom_highbd_obmc_variance128x64_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance64x128_c,
                        aom_highbd_obmc_variance64x128_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance64x64_c,
                        aom_highbd_obmc_variance64x64_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance64x32_c,
                        aom_highbd_obmc_variance64x32_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance32x64_c,
                        aom_highbd_obmc_variance32x64_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance32x32_c,
                        aom_highbd_obmc_variance32x32_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance32x16_c,
                        aom_highbd_obmc_variance32x16_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance16x32_c,
                        aom_highbd_obmc_variance16x32_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance16x16_c,
                        aom_highbd_obmc_variance16x16_sse4_1, 8),
  TestFuncs(aom_highbd_obmc_variance16x8_c, aom_highbd_obmc_variance16x8_sse4_1,
                        8),
  TestFuncs(aom_highbd_obmc_variance8x16_c, aom_highbd_obmc_variance8x16_sse4_1,
                        8),
  TestFuncs(aom_highbd_obmc_variance8x8_c, aom_highbd_obmc_variance8x8_sse4_1,
                        8),
  TestFuncs(aom_highbd_obmc_variance8x4_c, aom_highbd_obmc_variance8x4_sse4_1,
                        8),
  TestFuncs(aom_highbd_obmc_variance4x8_c, aom_highbd_obmc_variance4x8_sse4_1,
                        8),
  TestFuncs(aom_highbd_obmc_variance4x4_c, aom_highbd_obmc_variance4x4_sse4_1,
                        8),
  TestFuncs(aom_highbd_10_obmc_variance128x128_c,
                        aom_highbd_10_obmc_variance128x128_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance128x64_c,
                        aom_highbd_10_obmc_variance128x64_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance64x128_c,
                        aom_highbd_10_obmc_variance64x128_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance64x64_c,
                        aom_highbd_10_obmc_variance64x64_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance64x32_c,
                        aom_highbd_10_obmc_variance64x32_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance32x64_c,
                        aom_highbd_10_obmc_variance32x64_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance32x32_c,
                        aom_highbd_10_obmc_variance32x32_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance32x16_c,
                        aom_highbd_10_obmc_variance32x16_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance16x32_c,
                        aom_highbd_10_obmc_variance16x32_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance16x16_c,
                        aom_highbd_10_obmc_variance16x16_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance16x8_c,
                        aom_highbd_10_obmc_variance16x8_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance8x16_c,
                        aom_highbd_10_obmc_variance8x16_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance8x8_c,
                        aom_highbd_10_obmc_variance8x8_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance8x4_c,
                        aom_highbd_10_obmc_variance8x4_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance4x8_c,
                        aom_highbd_10_obmc_variance4x8_sse4_1, 10),
  TestFuncs(aom_highbd_10_obmc_variance4x4_c,
                        aom_highbd_10_obmc_variance4x4_sse4_1, 10),
  TestFuncs(aom_highbd_12_obmc_variance128x128_c,
                        aom_highbd_12_obmc_variance128x128_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance128x64_c,
                        aom_highbd_12_obmc_variance128x64_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance64x128_c,
                        aom_highbd_12_obmc_variance64x128_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance64x64_c,
                        aom_highbd_12_obmc_variance64x64_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance64x32_c,
                        aom_highbd_12_obmc_variance64x32_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance32x64_c,
                        aom_highbd_12_obmc_variance32x64_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance32x32_c,
                        aom_highbd_12_obmc_variance32x32_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance32x16_c,
                        aom_highbd_12_obmc_variance32x16_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance16x32_c,
                        aom_highbd_12_obmc_variance16x32_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance16x16_c,
                        aom_highbd_12_obmc_variance16x16_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance16x8_c,
                        aom_highbd_12_obmc_variance16x8_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance8x16_c,
                        aom_highbd_12_obmc_variance8x16_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance8x8_c,
                        aom_highbd_12_obmc_variance8x8_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance8x4_c,
                        aom_highbd_12_obmc_variance8x4_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance4x8_c,
                        aom_highbd_12_obmc_variance4x8_sse4_1, 12),
  TestFuncs(aom_highbd_12_obmc_variance4x4_c,
                        aom_highbd_12_obmc_variance4x4_sse4_1, 12)
};
*/
HBDObmcVarianceTestParam HIGHBD_TEST_PARAM[] = {
    HBDObmcVarianceTestParam(nullptr, nullptr, 8),
};

INSTANTIATE_TEST_CASE_P(OBMC, HBDObmcVarianceTest,
                        ::testing::ValuesIn(HIGHBD_TEST_PARAM));

}  // namespace
