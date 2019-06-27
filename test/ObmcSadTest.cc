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
 * @file ObmcSadTest.cc
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

typedef uint32_t (*EB_OBMC_SAD_FCN)(const uint8_t* pred, int pred_stride,
                                    const int32_t* wsrc, const int32_t* msk);

typedef std::tuple<EB_OBMC_SAD_FCN, EB_OBMC_SAD_FCN> ObmcSadTestParam;

class ObmcSadTest : public ::testing::TestWithParam<ObmcSadTestParam> {
  protected:
    void check_obmc_random() {
        EB_OBMC_SAD_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_SAD_FCN tst_fcn = TEST_GET_PARAM(1);
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
            const uint32_t ref_res = ref_fcn(pre, pre_stride, wsrc, mask);
            const uint32_t tst_res = tst_fcn(pre, pre_stride, wsrc, mask);
            ASSERT_EQ(tst_res, ref_res);
        }
        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }

    void check_obmc_exterme_values() {
        EB_OBMC_SAD_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_SAD_FCN tst_fcn = TEST_GET_PARAM(1);
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
            const uint32_t ref_res = ref_fcn(pre, pre_stride, wsrc, mask);
            const uint32_t tst_res = tst_fcn(pre, pre_stride, wsrc, mask);
            ASSERT_EQ(tst_res, ref_res);
        }
        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }
};

TEST_P(ObmcSadTest, RandomValues) {
    check_obmc_random();
}

TEST_P(ObmcSadTest, ExtremeValues) {
    check_obmc_exterme_values();
}

/* TODO: fill the target function and reference function in test vector
#if HAVE_SSE4_1
const ObmcSadTest::ParamType sse4_functions[] = {
  TestFuncs(aom_obmc_sad128x128_c, aom_obmc_sad128x128_sse4_1),
  TestFuncs(aom_obmc_sad128x64_c, aom_obmc_sad128x64_sse4_1),
  TestFuncs(aom_obmc_sad64x128_c, aom_obmc_sad64x128_sse4_1),
  TestFuncs(aom_obmc_sad64x64_c, aom_obmc_sad64x64_sse4_1),
  TestFuncs(aom_obmc_sad64x32_c, aom_obmc_sad64x32_sse4_1),
  TestFuncs(aom_obmc_sad32x64_c, aom_obmc_sad32x64_sse4_1),
  TestFuncs(aom_obmc_sad32x32_c, aom_obmc_sad32x32_sse4_1),
  TestFuncs(aom_obmc_sad32x16_c, aom_obmc_sad32x16_sse4_1),
  TestFuncs(aom_obmc_sad16x32_c, aom_obmc_sad16x32_sse4_1),
  TestFuncs(aom_obmc_sad16x16_c, aom_obmc_sad16x16_sse4_1),
  TestFuncs(aom_obmc_sad16x8_c, aom_obmc_sad16x8_sse4_1),
  TestFuncs(aom_obmc_sad8x16_c, aom_obmc_sad8x16_sse4_1),
  TestFuncs(aom_obmc_sad8x8_c, aom_obmc_sad8x8_sse4_1),
  TestFuncs(aom_obmc_sad8x4_c, aom_obmc_sad8x4_sse4_1),
  TestFuncs(aom_obmc_sad4x8_c, aom_obmc_sad4x8_sse4_1),
  TestFuncs(aom_obmc_sad4x4_c, aom_obmc_sad4x4_sse4_1)
};
*/
ObmcSadTestParam TEST_PARAM[] = {
    ObmcSadTestParam(nullptr, nullptr),
};

INSTANTIATE_TEST_CASE_P(OBMC, ObmcSadTest, ::testing::ValuesIn(TEST_PARAM));

////////////////////////////////////////////////////////////////////////////////
// High bit-depth
////////////////////////////////////////////////////////////////////////////////
class HBDObmcSadTest : public ::testing::TestWithParam<ObmcSadTestParam> {
  protected:
    void check_highbd_obmc_random() {
        EB_OBMC_SAD_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_SAD_FCN tst_fcn = TEST_GET_PARAM(1);
        ASSERT_NE(ref_fcn, nullptr) << "can not find reference function";
        ASSERT_NE(tst_fcn, nullptr) << "can not find test function";

        uint16_t* pre =
            (uint16_t*)aom_memalign(32, sizeof(uint16_t) * MAX_SB_SQUARE);
        int32_t* wsrc =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);
        int32_t* mask =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);

        SVTRandom rnd(12, false);
        SVTRandom rnd_stride(0, MAX_SB_SIZE + 1);
        SVTRandom rnd_mask(0, kMaskMax * kMaskMax + 1);
        for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
            for (int i = 0; i < MAX_SB_SQUARE; ++i) {
                pre[i] = rnd.random();
                wsrc[i] = rnd.random() * rnd_mask.random();
                mask[i] = rnd_mask.random();
            }

            const int pre_stride = rnd_stride.random();
            const uint32_t ref_res =
                ref_fcn((uint8_t*)pre, pre_stride, wsrc, mask);
            const uint32_t tst_res =
                tst_fcn((uint8_t*)pre, pre_stride, wsrc, mask);
            ASSERT_EQ(tst_res, ref_res);
        }
        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }

    void check_highbd_obmc_exterme_values() {
        EB_OBMC_SAD_FCN ref_fcn = TEST_GET_PARAM(0);
        EB_OBMC_SAD_FCN tst_fcn = TEST_GET_PARAM(1);
        ASSERT_NE(ref_fcn, nullptr) << "can not find reference function";
        ASSERT_NE(tst_fcn, nullptr) << "can not find test function";

        uint16_t* pre =
            (uint16_t*)aom_memalign(32, sizeof(uint16_t) * MAX_SB_SQUARE);
        int32_t* wsrc =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);
        int32_t* mask =
            (int32_t*)aom_memalign(32, sizeof(int32_t) * MAX_SB_SQUARE);

        for (int iter = 0; iter < MAX_SB_SIZE && !HasFatalFailure(); ++iter) {
            for (int i = 0; i < MAX_SB_SQUARE; ++i) {
                pre[i] = (1 << 12) - 1;
                wsrc[i] = ((1 << 12) - 1) * kMaskMax * kMaskMax;
                mask[i] = kMaskMax * kMaskMax;
            }

            const int pre_stride = iter;
            const uint32_t ref_res =
                ref_fcn((uint8_t*)pre, pre_stride, wsrc, mask);
            const uint32_t tst_res =
                tst_fcn((uint8_t*)pre, pre_stride, wsrc, mask);
            ASSERT_EQ(tst_res, ref_res);
        }

        aom_free(pre);
        aom_free(wsrc);
        aom_free(mask);
    }
};

TEST_P(HBDObmcSadTest, RandomValues) {
    check_highbd_obmc_random();
}

TEST_P(HBDObmcSadTest, ExtremeValues) {
    check_highbd_obmc_exterme_values();
}

/* TODO: fill the target function and reference function in test vector
#if HAVE_SSE4_1
ObmcSadHBDTest::ParamType sse4_functions_hbd[] = {
  TestFuncs(aom_highbd_obmc_sad128x128_c, aom_highbd_obmc_sad128x128_sse4_1),
  TestFuncs(aom_highbd_obmc_sad128x64_c, aom_highbd_obmc_sad128x64_sse4_1),
  TestFuncs(aom_highbd_obmc_sad64x128_c, aom_highbd_obmc_sad64x128_sse4_1),
  TestFuncs(aom_highbd_obmc_sad64x64_c, aom_highbd_obmc_sad64x64_sse4_1),
  TestFuncs(aom_highbd_obmc_sad64x32_c, aom_highbd_obmc_sad64x32_sse4_1),
  TestFuncs(aom_highbd_obmc_sad32x64_c, aom_highbd_obmc_sad32x64_sse4_1),
  TestFuncs(aom_highbd_obmc_sad32x32_c, aom_highbd_obmc_sad32x32_sse4_1),
  TestFuncs(aom_highbd_obmc_sad32x16_c, aom_highbd_obmc_sad32x16_sse4_1),
  TestFuncs(aom_highbd_obmc_sad16x32_c, aom_highbd_obmc_sad16x32_sse4_1),
  TestFuncs(aom_highbd_obmc_sad16x16_c, aom_highbd_obmc_sad16x16_sse4_1),
  TestFuncs(aom_highbd_obmc_sad16x8_c, aom_highbd_obmc_sad16x8_sse4_1),
  TestFuncs(aom_highbd_obmc_sad8x16_c, aom_highbd_obmc_sad8x16_sse4_1),
  TestFuncs(aom_highbd_obmc_sad8x8_c, aom_highbd_obmc_sad8x8_sse4_1),
  TestFuncs(aom_highbd_obmc_sad8x4_c, aom_highbd_obmc_sad8x4_sse4_1),
  TestFuncs(aom_highbd_obmc_sad4x8_c, aom_highbd_obmc_sad4x8_sse4_1),
  TestFuncs(aom_highbd_obmc_sad4x4_c, aom_highbd_obmc_sad4x4_sse4_1)
};
*/
ObmcSadTestParam HIGHBD_TEST_PARAM[] = {
    ObmcSadTestParam(nullptr, nullptr),
};

INSTANTIATE_TEST_CASE_P(OBMC, HBDObmcSadTest,
                        ::testing::ValuesIn(HIGHBD_TEST_PARAM));

}  // namespace
