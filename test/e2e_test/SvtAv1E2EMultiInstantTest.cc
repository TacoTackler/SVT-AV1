/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */
/******************************************************************************
 * @file SvtAv1E2EMultiInstantTest.cc
 *
 * @brief Impelmentation of multipule instant of encoder test at same time
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include <thread>
#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "RefDecoder.h"
#include "SvtAv1E2EFramework.h"
#include "../random.h"

using namespace svt_av1_e2e_test;
using namespace svt_av1_test_tool;

class EncTestCtxt {
  public:
    EncTestCtxt() = delete;
    EncTestCtxt(const TestVideoVector& vector, const uint32_t id,
                uint32_t max_ch) {
        video_src_ = SvtAv1E2ETestFramework::prepare_video_src(vector);
        memset(&av1enc_ctx_, 0, sizeof(av1enc_ctx_));
        recon_queue_ = nullptr;
        refer_dec_ = nullptr;
        obu_frame_header_size_ = 0;
        collect_ = nullptr;
        ref_compare_ = nullptr;
        channel_id_ = id;
        max_channel_num_ = max_ch;
        src_frame_count_ = 0;
        src_eos_ = false;
        enc_eos_ = false;
        recon_eos_ = false;
    }
    ~EncTestCtxt() {
        if (video_src_) {
            delete video_src_;
            video_src_ = nullptr;
        }
        if (recon_queue_) {
            delete recon_queue_;
            recon_queue_ = nullptr;
        }
        if (refer_dec_) {
            delete refer_dec_;
            refer_dec_ = nullptr;
        }
        if (collect_) {
            delete collect_;
            collect_ = nullptr;
        }
        if (ref_compare_) {
            delete ref_compare_;
            ref_compare_ = nullptr;
        }
    }

  public:
    void init_test() {
        EbErrorType return_error = EB_ErrorNone;

        // check for video source
        ASSERT_NE(video_src_, nullptr) << "video source create failed!";
        return_error = video_src_->open_source(0, 0);
        ASSERT_EQ(return_error, EB_ErrorNone)
            << "open_source return error:" << return_error;
        // Check input parameters
        uint32_t width = video_src_->get_width_with_padding();
        uint32_t height = video_src_->get_height_with_padding();
        uint32_t bit_depth = video_src_->get_bit_depth();
        ASSERT_GT(width, 0) << "Video vector width error.";
        ASSERT_GT(height, 0) << "Video vector height error.";
        ASSERT_TRUE(bit_depth == 10 || bit_depth == 8)
            << "Video vector bitDepth error.";
        src_frame_count_ = video_src_->get_frame_count();
        ASSERT_GT(src_frame_count_, 0) << "video frame count zero.";

        // Init handle
        return_error = eb_init_handle(
            &av1enc_ctx_.enc_handle, &av1enc_ctx_, &av1enc_ctx_.enc_params);
        ASSERT_EQ(return_error, EB_ErrorNone)
            << "eb_init_handle return error:" << return_error;
        ASSERT_NE(av1enc_ctx_.enc_handle, nullptr)
            << "eb_init_handle return null handle.";
        SvtAv1E2ETestFramework::setup_src_param(video_src_,
                                                av1enc_ctx_.enc_params);
        /** setup channel id for instant identify */
        av1enc_ctx_.enc_params.channel_id = channel_id_;
        av1enc_ctx_.enc_params.active_channel_count = max_channel_num_;

#if TILES
        EbBool has_tiles = (EbBool)(av1enc_ctx_.enc_params.tile_columns ||
                                    av1enc_ctx_.enc_params.tile_rows);
#else
        EbBool has_tiles = (EbBool)EB_FALSE;
#endif
        obu_frame_header_size_ =
            has_tiles ? OBU_FRAME_HEADER_SIZE + 1 : OBU_FRAME_HEADER_SIZE;

        // create recon sink before setup parameter of encoder
        VideoFrameParam param;
        memset(&param, 0, sizeof(param));
        param.format = video_src_->get_image_format();
        param.width = video_src_->get_width_with_padding();
        param.height = video_src_->get_height_with_padding();
        recon_queue_ = create_frame_queue(param);
        ASSERT_NE(recon_queue_, nullptr) << "can not create recon sink!!";
        recon_queue_->set_frame_count(src_frame_count_);
        av1enc_ctx_.enc_params.recon_enabled = 1;

        // set param and init encoder
        return_error = eb_svt_enc_set_parameter(av1enc_ctx_.enc_handle,
                                                &av1enc_ctx_.enc_params);
        ASSERT_EQ(return_error, EB_ErrorNone)
            << "eb_svt_enc_set_parameter return error:" << return_error;
        return_error = eb_init_encoder(av1enc_ctx_.enc_handle);
        ASSERT_EQ(return_error, EB_ErrorNone)
            << "eb_init_encoder return error:" << return_error;

        // Prepare buffer
        // Input Buffer
        av1enc_ctx_.input_picture_buffer = new EbBufferHeaderType;
        ASSERT_NE(av1enc_ctx_.input_picture_buffer, nullptr)
            << "Malloc memory for inputPictureBuffer failed.";
        av1enc_ctx_.input_picture_buffer->p_buffer = nullptr;
        av1enc_ctx_.input_picture_buffer->size = sizeof(EbBufferHeaderType);
        av1enc_ctx_.input_picture_buffer->p_app_private = nullptr;
        av1enc_ctx_.input_picture_buffer->pic_type = EB_AV1_INVALID_PICTURE;
        // Output buffer
        av1enc_ctx_.output_stream_buffer = new EbBufferHeaderType;
        ASSERT_NE(av1enc_ctx_.output_stream_buffer, nullptr)
            << "Malloc memory for outputStreamBuffer failed.";
        av1enc_ctx_.output_stream_buffer->p_buffer =
            new uint8_t[EB_OUTPUTSTREAMBUFFERSIZE_MACRO(width * height)];
        ASSERT_NE(av1enc_ctx_.output_stream_buffer->p_buffer, nullptr)
            << "Malloc memory for outputStreamBuffer->p_buffer failed.";
        av1enc_ctx_.output_stream_buffer->size = sizeof(EbBufferHeaderType);
        av1enc_ctx_.output_stream_buffer->n_alloc_len =
            EB_OUTPUTSTREAMBUFFERSIZE_MACRO(width * height);
        av1enc_ctx_.output_stream_buffer->p_app_private = nullptr;
        av1enc_ctx_.output_stream_buffer->pic_type = EB_AV1_INVALID_PICTURE;

        // create reference decoder
        refer_dec_ = create_reference_decoder();
        ASSERT_NE(refer_dec_, nullptr) << "can not create reference decoder!!";
    }

    void close_test() {
        EbErrorType return_error = EB_ErrorNone;

        // close encoder
        return_error = eb_deinit_encoder(av1enc_ctx_.enc_handle);
        ASSERT_EQ(return_error, EB_ErrorNone)
            << "eb_deinit_encoder return error:" << return_error;

        // Destruct the component
        return_error = eb_deinit_handle(av1enc_ctx_.enc_handle);
        ASSERT_EQ(return_error, EB_ErrorNone)
            << "eb_deinit_handle return error:" << return_error;
        av1enc_ctx_.enc_handle = nullptr;

        // Clear
        if (av1enc_ctx_.output_stream_buffer != nullptr) {
            if (av1enc_ctx_.output_stream_buffer->p_buffer != nullptr) {
                delete[] av1enc_ctx_.output_stream_buffer->p_buffer;
            }
            delete av1enc_ctx_.output_stream_buffer;
            av1enc_ctx_.output_stream_buffer = nullptr;
        }
        if (av1enc_ctx_.input_picture_buffer != nullptr) {
            delete av1enc_ctx_.input_picture_buffer;
            av1enc_ctx_.input_picture_buffer = nullptr;
        }

        ASSERT_NE(video_src_, nullptr);
        video_src_->close_source();
    }

    void run_test_one_step() {
        EbErrorType return_error = EB_ErrorNone;

        // read video source frame
        if (!src_eos_) {
            uint8_t* frame = (uint8_t*)video_src_->get_next_frame();
            if (frame) {
                // Fill in Buffers Header control data
                av1enc_ctx_.input_picture_buffer->p_buffer = frame;
                av1enc_ctx_.input_picture_buffer->n_filled_len =
                    video_src_->get_frame_size();
                av1enc_ctx_.input_picture_buffer->flags = 0;
                av1enc_ctx_.input_picture_buffer->p_app_private = nullptr;
                av1enc_ctx_.input_picture_buffer->pts =
                    video_src_->get_frame_index();
                av1enc_ctx_.input_picture_buffer->pic_type =
                    EB_AV1_INVALID_PICTURE;
                // Send the picture
                EXPECT_EQ(EB_ErrorNone,
                          return_error = eb_svt_enc_send_picture(
                              av1enc_ctx_.enc_handle,
                              av1enc_ctx_.input_picture_buffer))
                    << "eb_svt_enc_send_picture error at: "
                    << av1enc_ctx_.input_picture_buffer->pts;
            } else {
                src_eos_ = true;
                EbBufferHeaderType headerPtrLast;
                headerPtrLast.n_alloc_len = 0;
                headerPtrLast.n_filled_len = 0;
                headerPtrLast.n_tick_count = 0;
                headerPtrLast.p_app_private = nullptr;
                headerPtrLast.flags = EB_BUFFERFLAG_EOS;
                headerPtrLast.p_buffer = nullptr;
                headerPtrLast.pic_type = EB_AV1_INVALID_PICTURE;
                av1enc_ctx_.input_picture_buffer->flags = EB_BUFFERFLAG_EOS;
                EXPECT_EQ(EB_ErrorNone,
                          return_error = eb_svt_enc_send_picture(
                              av1enc_ctx_.enc_handle, &headerPtrLast))
                    << "eb_svt_enc_send_picture EOS error";
            }
        }

        // get recon
        if (recon_queue_ && !recon_eos_) {
            SvtAv1E2ETestFramework::get_recon_frame(
                av1enc_ctx_, recon_queue_, recon_eos_);
        }

        // get encoder output can check with recon
        if (!enc_eos_) {
            do {
                // non-blocking call
                EbBufferHeaderType* enc_out = nullptr;
                int pic_send_done = (src_eos_ && recon_eos_) ? 1 : 0;
                return_error = eb_svt_get_packet(
                    av1enc_ctx_.enc_handle, &enc_out, pic_send_done);
                ASSERT_NE(return_error, EB_ErrorMax)
                    << "Error while encoding, code:" << enc_out->flags;

                // process the output buffer
                if (return_error != EB_NoErrorEmptyQueue && enc_out) {
                    if (enc_out->flags & EB_BUFFERFLAG_SHOW_EXT) {
                        uint32_t first_part_size = enc_out->n_filled_len -
                                                   obu_frame_header_size_ -
                                                   TD_SIZE;
                        decode_compress_data(enc_out->p_buffer,
                                             first_part_size);
                        decode_compress_data(
                            enc_out->p_buffer + first_part_size,
                            obu_frame_header_size_ + TD_SIZE);
                    } else {
                        decode_compress_data(enc_out->p_buffer,
                                             enc_out->n_filled_len);
                    }
                    // encoder eos
                    if (enc_out->flags & EB_BUFFERFLAG_EOS) {
                        enc_eos_ = true;
                        printf("Encoder EOS\n");
                        break;
                    }
                } else {
                    if (return_error != EB_NoErrorEmptyQueue) {
                        enc_eos_ = true;
                        GTEST_FAIL() << "decoder return: " << return_error;
                    }
                    break;
                }

                // Release the output buffer
                if (enc_out != nullptr) {
                    eb_svt_release_out_buffer(&enc_out);
                }
            } while (src_eos_);
        }

        // flush un-checked frames in the end
        if (is_complete() && ref_compare_) {
            ref_compare_->flush_video();
            printf("test instance[%u] is complete.\n", channel_id_);
        }
    }

    bool is_complete() {
        return src_eos_ && recon_eos_ && enc_eos_;
    }

  private:
    void decode_compress_data(const uint8_t* data, const uint32_t size) {
        ASSERT_NE(data, nullptr);
        ASSERT_GT(size, 0);

        // input the compressed data into decoder
        ASSERT_EQ(refer_dec_->decode(data, size), RefDecoder::REF_CODEC_OK);

        VideoFrame ref_frame;
        memset(&ref_frame, 0, sizeof(ref_frame));
        while (refer_dec_->get_frame(ref_frame) == RefDecoder::REF_CODEC_OK) {
            if (recon_queue_) {
                // compare tools
                if (ref_compare_ == nullptr) {
                    ref_compare_ =
                        create_ref_compare_queue(ref_frame, recon_queue_);
                    ASSERT_NE(ref_compare_, nullptr);
                }
                // Compare ref decode output with recon output.
                ASSERT_TRUE(ref_compare_->compare_video(ref_frame))
                    << "image compare failed on " << ref_frame.timestamp;
            }
        }
    }

  protected:
    uint32_t channel_id_;           /**< channel id for encoder instance */
    uint32_t max_channel_num_;      /**< maximum channel available */
    VideoSource* video_src_;        /**< video source context */
    SvtAv1Context av1enc_ctx_;      /**< AV1 encoder context */
    FrameQueue* recon_queue_;       /**< reconstruction frame collection */
    RefDecoder* refer_dec_;         /**< reference decoder context */
    uint8_t obu_frame_header_size_; /**< size of obu frame header */
    PerformanceCollect* collect_;   /**< performance and time collection*/
    ICompareQueue* ref_compare_; /**< sink of reference to compare with recon*/
    uint32_t src_frame_count_;   /**< frame counter of source video*/
    bool src_eos_;               /**< end flag of video source */
    bool recon_eos_;             /**< end flag of recon frames */
    bool enc_eos_;               /**< end flag of encoder output */
};

class SvtAv1E2EMultiInstSerialTest
    : public ::testing::TestWithParam<MultiInstVector> {
  public:
    SvtAv1E2EMultiInstSerialTest() : rand_(0, 99) {
        inst_vec_.clear();
        num_inst_ = std::get<1>(GetParam());
    }
    virtual ~SvtAv1E2EMultiInstSerialTest() {
    }

  protected:
    void SetUp() override {
        for (uint32_t i = 0; i < num_inst_; ++i) {
            EncTestCtxt* new_inst =
                new EncTestCtxt(std::get<0>(GetParam()), i, num_inst_);
            ASSERT_NE(new_inst, nullptr) << "create new test instance failed.";
            new_inst->init_test();
            inst_vec_.push_back(new_inst);
        }
    }
    void TearDown() override {
        while (inst_vec_.size()) {
            EncTestCtxt* inst = inst_vec_.back();
            inst->close_test();
            delete inst;
            inst_vec_.pop_back();
        }
    }
    bool is_all_complete() {
        for (EncTestCtxt* inst : inst_vec_) {
            if (inst && !inst->is_complete()) {
                return false;
            }
        }
        return true;
    }
    void run_serialize_check() {
        while (!is_all_complete()) {
            for (EncTestCtxt* inst : inst_vec_) {
                ASSERT_NE(inst, nullptr);
                if (!inst->is_complete())
                    inst->run_test_one_step();
            }
        }
    }
    void run_random_serialize_check() {
        while (!is_all_complete()) {
            // random get test instance invector
            EncTestCtxt* inst = inst_vec_.at(rand_.random() % num_inst_);
            ASSERT_NE(inst, nullptr);
            if (!inst->is_complete())
                inst->run_test_one_step();
        }
    }

  private:
    uint32_t num_inst_;                  /**< number of instances */
    std::vector<EncTestCtxt*> inst_vec_; /**< vector of test instances */
    SVTRandom rand_;                     /**< random generator */
};

TEST_P(SvtAv1E2EMultiInstSerialTest, run_serialize_check) {
    run_serialize_check();
}

TEST_P(SvtAv1E2EMultiInstSerialTest, run_random_serialize_check) {
    run_random_serialize_check();
}

INSTANTIATE_TEST_CASE_P(SVT_AV1, SvtAv1E2EMultiInstSerialTest,
                        ::testing::ValuesIn(multi_inst_vectors));

class SvtAv1E2EMultiInstParallelTest
    : public ::testing::TestWithParam<MultiInstVector> {
  public:
    typedef struct InstanceBind {
        std::thread* run_thrd;
        EncTestCtxt* test_ctxt;
    } InstanceBind;

  public:
    SvtAv1E2EMultiInstParallelTest() {
        inst_vec_.clear();
        num_inst_ = std::get<1>(GetParam());
    }
    virtual ~SvtAv1E2EMultiInstParallelTest() {
    }

  protected:
    void SetUp() override {
        for (uint32_t i = 0; i < num_inst_; ++i) {
            EncTestCtxt* new_inst =
                new EncTestCtxt(std::get<0>(GetParam()), i, num_inst_);
            ASSERT_NE(new_inst, nullptr) << "create new test instance failed.";
            new_inst->init_test();
            std::thread* new_thrd =
                new std::thread(run_instance_check, new_inst);
            InstanceBind inst_bind = {new_thrd, new_inst};
            inst_vec_.push_back(inst_bind);
        }
    }
    void TearDown() override {
        while (inst_vec_.size()) {
            InstanceBind inst_bind = inst_vec_.back();
            inst_bind.test_ctxt->close_test();
            inst_bind.run_thrd->join();
            delete inst_bind.run_thrd;
            delete inst_bind.test_ctxt;
            inst_vec_.pop_back();
        }
    }
    bool is_all_complete() {
        for (InstanceBind inst : inst_vec_) {
            if (inst.test_ctxt && !inst.test_ctxt->is_complete()) {
                return false;
            }
        }
        return true;
    }
    void run_parallel_check() {
        while (!is_all_complete()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
    static void run_instance_check(EncTestCtxt* inst) {
        ASSERT_NE(inst, nullptr);
        while (!inst->is_complete()) {
            inst->run_test_one_step();
        }
    }

  private:
    uint32_t num_inst_;                  /**< number of instances */
    std::vector<InstanceBind> inst_vec_; /**< vector of test instances */
};

TEST_P(SvtAv1E2EMultiInstParallelTest, run_parallel_check) {
    run_parallel_check();
}

INSTANTIATE_TEST_CASE_P(SVT_AV1, SvtAv1E2EMultiInstParallelTest,
                        ::testing::ValuesIn(multi_inst_vectors));
