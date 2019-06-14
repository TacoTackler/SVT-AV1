/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SvtAv1E2ETest.cc
 *
 * @brief Impelmentation of SVT-AV1 encoder E2E test
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "EbSvtAv1Enc.h"
#include "gtest/gtest.h"
#include "SvtAv1E2EFramework.h"
#include "ConfigReader.h"

using namespace svt_av1_e2e_test;
using namespace svt_av1_e2e_test_vector;

/**
 * @brief SVT-AV1 encoder simple E2E test
 *
 * Test strategy:
 * Setup SVT-AV1 encoder with default parameter, and encode the input YUV data
 * frames.
 *
 * Expected result:
 * No error is reported in encoding progress. The output compressed data
 * is complete.
 *
 * Test coverage:
 * All test vectors
 */
class SvtAv1E2ESimpleTest : public SvtAv1E2ETestFramework {};

TEST_P(SvtAv1E2ESimpleTest, run_smoking_test) {
    run_encode_process();
}

INSTANTIATE_TEST_CASE_P(
    SVT_AV1, SvtAv1E2ESimpleTest,
    ::testing::ValuesIn(generate_vector_from_config("video_src.cfg")));

/**
 * @brief SVT-AV1 encoder simple E2E test with save compressed data in file
 *
 * Test strategy:
 * Setup SVT-AV1 encoder with default parameter, and encode the input YUV data
 * frames. Save the compressed data into IVF file.
 *
 * Expected result:
 * No error is reported in encoding progress. The output compressed data
 * is saved into IVF file.
 *
 * Test coverage:
 * Smoking test vectors
 */
class SvtAv1E2ESimpleFileTest : public SvtAv1E2ETestFramework {
  protected:
    /** initialization for test */
    void init_test() override {
        output_file_ = new IvfFile("output.av1");
        SvtAv1E2ETestFramework::init_test();
    }
};

TEST_P(SvtAv1E2ESimpleFileTest, run_smoking_with_output_test) {
    run_encode_process();
}

INSTANTIATE_TEST_CASE_P(
    SVT_AV1, SvtAv1E2ESimpleFileTest,
    ::testing::ValuesIn(generate_vector_from_config("smoking_test.cfg")));

/**
 * @brief SVT-AV1 encoder E2E test with comparing the reconstructed frame with
 * output frame from decoder buffer list
 *
 * Test strategy:
 * Setup SVT-AV1 encoder with default parameter, and encode the input YUV data
 * frames. Collect the reconstructed frames and compared them with reference
 * decoder output.
 *
 * Expected result:
 * No error is reported in encoding progress. The reconstructed frame
 * data is same as the output frame from reference decoder.
 *
 * Test coverage:
 * All test vectors
 */
class SvtAv1E2EConformanceTest : public SvtAv1E2ETestFramework {
  protected:
    /** initialization for test */
    void init_test() override {
        // create recon sink before setup parameter of encoder
        VideoFrameParam param;
        memset(&param, 0, sizeof(param));
        param.format = video_src_->get_image_format();
        param.width = video_src_->get_width_with_padding();
        param.height = video_src_->get_height_with_padding();
        recon_queue_ = create_frame_queue(param);
        ASSERT_NE(recon_queue_, nullptr) << "can not create recon sink!!";
        if (recon_queue_)
            av1enc_ctx_.enc_params.recon_enabled = 1;

        // create reference decoder
        refer_dec_ = create_reference_decoder();
        ASSERT_NE(refer_dec_, nullptr) << "can not create reference decoder!!";

        collect_ = new PerformanceCollect(typeid(this).name());

        SvtAv1E2ETestFramework::init_test();
    }
};

TEST_P(SvtAv1E2EConformanceTest, run_conformance_test) {
    run_encode_process();
}

INSTANTIATE_TEST_CASE_P(
    SVT_AV1, SvtAv1E2EConformanceTest,
    ::testing::ValuesIn(generate_vector_from_config("conformance_test.cfg")));

#define COPY_PARAM_VALUE(param_name)                                \
    {                                                               \
        bool can_copy = false;                                      \
        if (sizeof(config->param_name) == 1 &&                      \
            ((uint8_t)config->param_name) != 0xFF)                  \
            can_copy = true;                                        \
        else if (sizeof(config->param_name) == 2 &&                 \
                 ((uint16_t)config->param_name) != 0xFFFF)          \
            can_copy = true;                                        \
        else if (sizeof(config->param_name) == 4 &&                 \
                 ((uint32_t)config->param_name) != 0xFFFFFFFF)      \
            can_copy = true;                                        \
        else if (sizeof(config->param_name) == 8 &&                 \
                 config->param_name != 0xFFFFFFFFFFFFFFFF)          \
            can_copy = true;                                        \
        if (can_copy)                                               \
            av1enc_ctx_.enc_params.param_name = config->param_name; \
    }

class SvtAv1E2EConfigurableTest : public SvtAv1E2EConformanceTest {
  public:
    SvtAv1E2EConfigurableTest() {
        std::string config_path = VideoFileSource::get_vector_dir();
        config_path += "/e2e_test.cfg";
        reader_ = new ConfigReader(config_path);
        test_case_index_ = 0;
    }
    ~SvtAv1E2EConfigurableTest() {
    }

    /** initialization for test */
    void init_test() override {
        ASSERT_NE(reader_, nullptr);
        const ConfigReader::TestConfig* config =
            reader_->get_config(test_case_index_);
        if (config == nullptr)
            return;
        printf("--------- Configuration %d ----------\n", test_case_index_);
        config->print_param();
        printf("---------      End %d      ----------\n", test_case_index_);

        copy_param_setting(&config->param_config);

        /** create recon frame queue before setup parameter of encoder */
        VideoFrameParam param;
        memset(&param, 0, sizeof(param));
        param.format = video_src_->get_image_format();
        param.width = video_src_->get_width_with_padding();
        param.height = video_src_->get_height_with_padding();
        recon_queue_ = create_frame_queue(param);
        ASSERT_NE(recon_queue_, nullptr) << "can not create recon queue!!";
        if (recon_queue_)
            av1enc_ctx_.enc_params.recon_enabled = 1;

        /** create reference decoder*/
        refer_dec_ = create_reference_decoder(config->use_parser);
        ASSERT_NE(refer_dec_, nullptr) << "can not create reference decoder!!";
        refer_dec_->set_resolution(video_src_->get_width_with_padding(),
                                   video_src_->get_height_with_padding());

        SvtAv1E2ETestFramework::init_test();
    }
    void close_test() override {
        SvtAv1E2ETestFramework::close_test();
        test_case_index_++;
    }

  protected:
    void copy_param_setting(const EbSvtAv1EncConfiguration* config) {
        COPY_PARAM_VALUE(enc_mode);
        COPY_PARAM_VALUE(intra_period_length);
        COPY_PARAM_VALUE(intra_refresh_type);
        COPY_PARAM_VALUE(hierarchical_levels);
        COPY_PARAM_VALUE(pred_structure);
        COPY_PARAM_VALUE(base_layer_switch_mode);
        COPY_PARAM_VALUE(source_width);
        COPY_PARAM_VALUE(source_height);
        COPY_PARAM_VALUE(frame_rate);
        COPY_PARAM_VALUE(frame_rate_numerator);
        COPY_PARAM_VALUE(frame_rate_denominator);
        COPY_PARAM_VALUE(encoder_bit_depth);
        COPY_PARAM_VALUE(encoder_color_format);
        COPY_PARAM_VALUE(compressed_ten_bit_format);
        COPY_PARAM_VALUE(frames_to_be_encoded);
        COPY_PARAM_VALUE(improve_sharpness);
        COPY_PARAM_VALUE(sb_sz);
        COPY_PARAM_VALUE(super_block_size);
        COPY_PARAM_VALUE(partition_depth);
        COPY_PARAM_VALUE(qp);
        COPY_PARAM_VALUE(use_qp_file);
        COPY_PARAM_VALUE(enable_qp_scaling_flag);
        COPY_PARAM_VALUE(disable_dlf_flag);
        COPY_PARAM_VALUE(enable_denoise_flag);
        COPY_PARAM_VALUE(film_grain_denoise_strength);
        COPY_PARAM_VALUE(enable_warped_motion);
        COPY_PARAM_VALUE(use_default_me_hme);
        COPY_PARAM_VALUE(enable_hme_flag);
        COPY_PARAM_VALUE(ext_block_flag);
        COPY_PARAM_VALUE(in_loop_me_flag);
        COPY_PARAM_VALUE(search_area_width);
        COPY_PARAM_VALUE(search_area_height);
        COPY_PARAM_VALUE(constrained_intra);
        COPY_PARAM_VALUE(rate_control_mode);
        COPY_PARAM_VALUE(scene_change_detection);
        COPY_PARAM_VALUE(look_ahead_distance);
        COPY_PARAM_VALUE(target_bit_rate);
        COPY_PARAM_VALUE(max_qp_allowed);
        COPY_PARAM_VALUE(min_qp_allowed);
        COPY_PARAM_VALUE(high_dynamic_range_input);
        COPY_PARAM_VALUE(profile);
        COPY_PARAM_VALUE(tier);
        COPY_PARAM_VALUE(level);
        COPY_PARAM_VALUE(asm_type);
        COPY_PARAM_VALUE(channel_id);
        COPY_PARAM_VALUE(active_channel_count);
        COPY_PARAM_VALUE(speed_control_flag);
        COPY_PARAM_VALUE(injector_frame_rate);
        COPY_PARAM_VALUE(logical_processors);
        COPY_PARAM_VALUE(target_socket);
        COPY_PARAM_VALUE(recon_enabled);
#if TILES
        COPY_PARAM_VALUE(tile_columns);
        COPY_PARAM_VALUE(tile_rows);
#endif
    }

    void run_all_configurable_test() {
        for (test_case_index_ = 0;
             test_case_index_ < reader_->get_config_count();
             test_case_index_++) {
            SvtAv1E2EConformanceTest::SetUp();
            run_encode_process();
            SvtAv1E2EConformanceTest::TearDown();
        }
    }
    void SetUp() override {
        /* skip SvtAv1E2EConformanceTest::SetUp() */
    }
    void TearDown() override {
        /* skip SvtAv1E2EConformanceTest::TearDown() */
    }

  protected:
    ConfigReader* reader_;     /**< configuration reader */
    uint32_t test_case_index_; /**< index of test case is running */
};

TEST_P(SvtAv1E2EConfigurableTest, run_conformance_test) {
    run_all_configurable_test();
}

INSTANTIATE_TEST_CASE_P(
    SVT_AV1, SvtAv1E2EConfigurableTest,
    ::testing::ValuesIn(generate_vector_from_config("conformance_test.cfg")));
