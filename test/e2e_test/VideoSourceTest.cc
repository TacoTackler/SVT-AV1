/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file VideoSourceTest.cc
 *
 * @brief
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "EbSvtAv1Enc.h"
#include "Y4mVideoSource.h"
#include "YuvVideoSource.h"
#include "DummyVideoSource.h"
#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "RefDecoder.h"
#include "SvtAv1E2EFramework.h"
#include "CompareTools.h"
#include "ConfigEncoder.h"

using namespace svt_av1_video_source;
using namespace svt_av1_e2e_test_vector;

namespace video_src_test {

using VideoSrcTestParam =
    std::tuple<const TestVideoVector, const TestVideoVector>;

class VideoSrcTest : public ::testing::TestWithParam<VideoSrcTestParam> {
  protected:
    VideoSrcTest()
        : vector1_(std::get<0>(GetParam())),
          vector2_(std::get<1>(GetParam())),
          test_src1(nullptr),
          test_src2(nullptr) {
    }
    virtual ~VideoSrcTest() {
    }

    void SetUp() override {
        test_src1 = prepare_video_src(vector1_);
        ASSERT_NE(test_src1, nullptr) << "video source 1 created failed!";
        if (std::get<0>(vector2_).size() > 0) {
            test_src2 = prepare_video_src(vector2_);
            ASSERT_NE(test_src2, nullptr) << "video source 2 created failed!";
        }
    }
    void TearDown() override {
        if (test_src1) {
            delete test_src1;
            test_src1 = nullptr;
        }
        if (test_src2) {
            delete test_src2;
            test_src2 = nullptr;
        }
    }

    void run_open_close_test(const int times) {
        EbErrorType err;
        for (int i = 0; i < times; i++) {
            // test for open source
            err = test_src1->open_source(std::get<7>(vector1_),
                                         std::get<8>(vector1_));
            ASSERT_EQ(err, EB_ErrorNone)
                << "video source 1 open_source return error:" << err;
            if (test_src2) {
                err = test_src2->open_source(std::get<7>(vector2_),
                                             std::get<8>(vector2_));
                ASSERT_EQ(err, EB_ErrorNone)
                    << "video source 1 open_source return error:" << err;
            }

            if (test_src2) {
                // compare two info of video source
                ASSERT_EQ(test_src1->get_width_with_padding(),
                          test_src2->get_width_with_padding());
                ASSERT_EQ(test_src1->get_height_with_padding(),
                          test_src2->get_height_with_padding());
                ASSERT_EQ(test_src1->get_bit_depth(),
                          test_src2->get_bit_depth());
                ASSERT_EQ(test_src1->get_compressed_10bit_mode(),
                          test_src2->get_compressed_10bit_mode());
                ASSERT_EQ(test_src1->get_frame_count(),
                          test_src2->get_frame_count());
                ASSERT_EQ(test_src1->get_image_format(),
                          test_src2->get_image_format());
                ASSERT_EQ(test_src1->get_frame_size(),
                          test_src2->get_frame_size());
                ASSERT_EQ(test_src1->get_frame_index(),
                          test_src2->get_frame_index());
            } else {
                // compare info with the value in vector
                uint32_t width = std::get<3>(vector1_);
                width = (width % 16 != 0) ? ((width >> 4) + 1) << 4 : width;
                ASSERT_EQ(test_src1->get_width_with_padding(), width);
                uint32_t height = std::get<4>(vector1_);
                height = (height % 16 != 0) ? ((height >> 4) + 1) << 4 : height;
                ASSERT_EQ(test_src1->get_height_with_padding(), height);
                ASSERT_EQ(test_src1->get_bit_depth(), std::get<5>(vector1_));
                ASSERT_EQ(test_src1->get_compressed_10bit_mode(),
                          std::get<6>(vector1_));
            }

            // test for close source
            test_src1->close_source();
            if (test_src2) {
                test_src2->close_source();
            }
        }
    }

    void run_read_frame_test() {
        EbErrorType err;
        // test for open source
        err = test_src1->open_source(std::get<7>(vector1_),
                                     std::get<8>(vector1_));
        ASSERT_EQ(err, EB_ErrorNone)
            << "video source 1 open_source return error:" << err;
        if (test_src2) {
            err = test_src2->open_source(std::get<7>(vector2_),
                                         std::get<8>(vector2_));
            ASSERT_EQ(err, EB_ErrorNone)
                << "video source 1 open_source return error:" << err;
        }

        EbSvtIOFormat* frame1 = nullptr;
        EbSvtIOFormat* frame2 = nullptr;
        uint32_t frame_count = 0;
        do {
            frame1 = test_src1->get_next_frame();
            if (test_src2) {
                frame2 = test_src2->get_next_frame();
                if (frame1 == nullptr) {
                    ASSERT_EQ(frame2, nullptr)
                        << "frame count mismatch in video source 1 and 2!";
                } else {
                    ASSERT_NE(frame2, nullptr)
                        << "can not get next frame in video source 2!";
                    compare_frames(frame1, frame2);
                }
            } else {
            }
            if (frame1 != nullptr) {
                frame_count++;
            }
        } while (frame1 && !HasFailure());
        ASSERT_EQ(frame_count, std::get<8>(vector1_))
            << "frame count of video source 1 is mismatched!";

        // test for close source
        test_src1->close_source();
        if (test_src2) {
            test_src2->close_source();
        }
    }

    static VideoSource* prepare_video_src(const TestVideoVector& vector) {
        VideoSource* video_src = nullptr;
        switch (std::get<1>(vector)) {
        case YUV_VIDEO_FILE:
            video_src = new YuvVideoSource(std::get<0>(vector),
                                           std::get<2>(vector),
                                           std::get<3>(vector),
                                           std::get<4>(vector),
                                           (uint8_t)std::get<5>(vector));
            break;
        case Y4M_VIDEO_FILE:
            video_src = new Y4MVideoSource(std::get<0>(vector),
                                           std::get<2>(vector),
                                           std::get<3>(vector),
                                           std::get<4>(vector),
                                           (uint8_t)std::get<5>(vector),
                                           std::get<6>(vector));
            break;
        case DUMMY_SOURCE:
            video_src = new DummyVideoSource(std::get<2>(vector),
                                             std::get<3>(vector),
                                             std::get<4>(vector),
                                             (uint8_t)std::get<5>(vector),
                                             std::get<6>(vector));
            break;
        default: assert(0); break;
        }
        return video_src;
    }

    void compare_frames(const EbSvtIOFormat* frame1,
                        const EbSvtIOFormat* frame2) {
        ASSERT_EQ(frame1->width, frame2->width);
        ASSERT_EQ(frame1->height, frame2->height);
        ASSERT_EQ(frame1->origin_x, frame2->origin_x);
        ASSERT_EQ(frame1->origin_y, frame2->origin_y);
        ASSERT_EQ(frame1->y_stride, frame2->y_stride);
        ASSERT_EQ(frame1->cb_stride, frame2->cb_stride);
        ASSERT_EQ(frame1->cr_stride, frame2->cr_stride);

        uint32_t width_downsize = 0;
        uint32_t height_downsize = 0;
        switch (std::get<2>(vector1_)) {
        case IMG_FMT_420P10_PACKED:
        case IMG_FMT_420: {
            width_downsize = 1;
            height_downsize = 1;
        } break;
        case IMG_FMT_422P10_PACKED:
        case IMG_FMT_422: {
            width_downsize = 1;
        } break;
        case IMG_FMT_444P10_PACKED:
        case IMG_FMT_444: break;
        default: assert(0); break;
        }

        uint32_t width = frame1->width;
        uint32_t height = frame1->height;
        uint32_t bytes_per_sample = std::get<5>(vector1_) > 8 ? 2 : 1;

        // luma compare
        uint8_t* pos1 = frame1->luma;
        uint8_t* pos2 = frame2->luma;
        uint32_t stride = frame1->y_stride * bytes_per_sample;
        for (size_t i = 0; i < height; i++) {
            for (size_t j = 0; j < stride; j++) {
                ASSERT_EQ(pos1[j], pos2[j])
                    << "luma compare failed at: [" << i << ":" << j << "]";
            }
            pos1 += stride;
            pos2 += stride;
        }

        // cb compare
        pos1 = frame1->cb;
        pos2 = frame2->cb;
        stride = frame1->cb_stride * bytes_per_sample;
        for (size_t i = 0; i < (height >> height_downsize); i++) {
            for (size_t j = 0; j < stride; j++) {
                ASSERT_EQ(pos1[j], pos2[j])
                    << "cb compare failed at: [" << i << ":" << j << "]";
            }
            pos1 += stride;
            pos2 += stride;
        }

        // cr compare
        pos1 = frame1->cr;
        pos2 = frame2->cr;
        stride = frame1->cr_stride * bytes_per_sample;
        for (size_t i = 0; i < (height >> height_downsize); i++) {
            for (size_t j = 0; j < stride; j++) {
                ASSERT_EQ(pos1[j], pos2[j])
                    << "cr compare failed at: [" << i << ":" << j << "]";
            }
            pos1 += stride;
            pos2 += stride;
        }
    }

  private:
    const TestVideoVector& vector1_;
    const TestVideoVector& vector2_;
    VideoSource* test_src1;
    VideoSource* test_src2;
};  // namespace video_src_test

TEST_P(VideoSrcTest, RunOpenCloseTest) {
    run_open_close_test(100);
}

TEST_P(VideoSrcTest, RunReadFrameTest) {
    run_read_frame_test();
}

const std::vector<TestVideoVector> basic_src_vectors = {
    std::make_tuple("colorbar_480p_8_420", DUMMY_SOURCE, IMG_FMT_420, 640, 480,
                    8, 0, 0, 100),
    std::make_tuple("colorbar_4k_10_420", DUMMY_SOURCE, IMG_FMT_420, 4096, 2160,
                    10, 0, 0, 60),
    std::make_tuple("kirland_640_480_30.yuv", YUV_VIDEO_FILE, IMG_FMT_420, 640,
                    480, 8, 0, 0, 60),
    std::make_tuple("park_joy_90p_8_420.y4m", Y4M_VIDEO_FILE, IMG_FMT_420, 160,
                    90, 8, 0, 0, 0),
    std::make_tuple("park_joy_90p_10_420.y4m", Y4M_VIDEO_FILE, IMG_FMT_420, 160,
                    90, 10, 0, 0, 0),
};

const std::vector<VideoSrcTestParam> generate_single_src(
    const std::vector<TestVideoVector>& vectors) {
    std::vector<VideoSrcTestParam> params_list;
    for (TestVideoVector test_vector : vectors) {
        params_list.push_back(VideoSrcTestParam(
            test_vector,
            std::make_tuple(
                "", DUMMY_SOURCE, IMG_FMT_YV12, 0, 0, 0, false, 0, 0)));
    }
    return params_list;
}

const std::vector<VideoSrcTestParam> compare_frame_vectors = {
    // 8-bit YUV420
    std::make_tuple(std::make_tuple("colorbar_480p_8_420", DUMMY_SOURCE,
                                    IMG_FMT_420, 640, 480, 8, 0, 0, 100),
                    std::make_tuple("colorbar_480p_8_420.yuv", YUV_VIDEO_FILE,
                                    IMG_FMT_420, 640, 480, 8, 0, 0, 100)),
    std::make_tuple(std::make_tuple("colorbar_480p_8_420", DUMMY_SOURCE,
                                    IMG_FMT_420, 640, 480, 8, 0, 0, 100),
                    std::make_tuple("colorbar_480p_8_420.y4m", Y4M_VIDEO_FILE,
                                    IMG_FMT_420, 640, 480, 8, 0, 0, 100)),
    // 8-bit YUV422
    std::make_tuple(std::make_tuple("colorbar_480p_8_422", DUMMY_SOURCE,
                                    IMG_FMT_422, 640, 480, 8, 0, 0, 100),
                    std::make_tuple("colorbar_480p_8_422.yuv", YUV_VIDEO_FILE,
                                    IMG_FMT_422, 640, 480, 8, 0, 0, 100)),
    std::make_tuple(std::make_tuple("colorbar_480p_8_422", DUMMY_SOURCE,
                                    IMG_FMT_422, 640, 480, 8, 0, 0, 100),
                    std::make_tuple("colorbar_480p_8_422.y4m", Y4M_VIDEO_FILE,
                                    IMG_FMT_422, 640, 480, 8, 0, 0, 100)),
    // 8-bit YUV444
    std::make_tuple(std::make_tuple("colorbar_480p_8_444", DUMMY_SOURCE,
                                    IMG_FMT_444, 640, 480, 8, 0, 0, 100),
                    std::make_tuple("colorbar_480p_8_444.yuv", YUV_VIDEO_FILE,
                                    IMG_FMT_444, 640, 480, 8, 0, 0, 100)),
    std::make_tuple(std::make_tuple("colorbar_480p_8_444", DUMMY_SOURCE,
                                    IMG_FMT_444, 640, 480, 8, 0, 0, 100),
                    std::make_tuple("colorbar_480p_8_444.y4m", Y4M_VIDEO_FILE,
                                    IMG_FMT_444, 640, 480, 8, 0, 0, 100)),
    // 10-bit YUV420
    std::make_tuple(std::make_tuple("colorbar_480p_10_420", DUMMY_SOURCE,
                                    IMG_FMT_420, 640, 480, 10, 0, 0, 100),
                    std::make_tuple("colorbar_480p_10_420.yuv", YUV_VIDEO_FILE,
                                    IMG_FMT_420, 640, 480, 10, 0, 0, 100)),
    std::make_tuple(std::make_tuple("colorbar_480p_10_420", DUMMY_SOURCE,
                                    IMG_FMT_420, 640, 480, 10, 0, 0, 100),
                    std::make_tuple("colorbar_480p_10_420.y4m", Y4M_VIDEO_FILE,
                                    IMG_FMT_420, 640, 480, 10, 0, 0, 100)),
    // 10-bit YUV422
    std::make_tuple(std::make_tuple("colorbar_480p_10_422", DUMMY_SOURCE,
                                    IMG_FMT_422, 640, 480, 10, 0, 0, 100),
                    std::make_tuple("colorbar_480p_10_422.yuv", YUV_VIDEO_FILE,
                                    IMG_FMT_422, 640, 480, 10, 0, 0, 100)),
    std::make_tuple(std::make_tuple("colorbar_480p_10_422", DUMMY_SOURCE,
                                    IMG_FMT_422, 640, 480, 10, 0, 0, 100),
                    std::make_tuple("colorbar_480p_10_422.y4m", Y4M_VIDEO_FILE,
                                    IMG_FMT_422, 640, 480, 10, 0, 0, 100)),
    // 10-bit YUV444
    std::make_tuple(std::make_tuple("colorbar_480p_10_444", DUMMY_SOURCE,
                                    IMG_FMT_444, 640, 480, 10, 0, 0, 100),
                    std::make_tuple("colorbar_480p_10_444.yuv", YUV_VIDEO_FILE,
                                    IMG_FMT_444, 640, 480, 10, 0, 0, 100)),
    std::make_tuple(std::make_tuple("colorbar_480p_10_444", DUMMY_SOURCE,
                                    IMG_FMT_444, 640, 480, 10, 0, 0, 100),
                    std::make_tuple("colorbar_480p_10_444.y4m", Y4M_VIDEO_FILE,
                                    IMG_FMT_444, 640, 480, 10, 0, 0, 100)),
};

INSTANTIATE_TEST_CASE_P(VSRC, VideoSrcTest,
                        ::testing::ValuesIn(compare_frame_vectors));

}  // namespace video_src_test
