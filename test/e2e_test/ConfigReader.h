/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file ConfigReader.h
 *
 * @brief Defines a tool of reading configuration from file to setup encoder.
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#ifndef _CONFIG_READER_H_
#define _CONFIG_READER_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include "EbSvtAv1Enc.h"

using std::string;
using std::vector;

#define MAX_CONFIG_STR_LEN (1024)

#define DEFINE_SET_PARAM(param_name)                                  \
    static void set_##param_name(const string value, TestConfig* p) { \
        p->param_config.param_name = std::stoul(value);               \
    }

#define MAP_GEN(param_name) \
    { #param_name, ConfigReader::set_##param_name }

class ConfigReader {
  public:
    typedef struct TestConfig {
        string name;        /**< name of test case in string */
        string description; /**< description of test case in string */
        string test_vector; /**< test vector for this case */
        bool use_parser;    /**< flag to enable/disable parser */
        bool is_death;      /**< flag to identify if this test is a death test*/
        string result;      /**< expected test result */
        string config_line; /**< configuration in a line */
        vector<string> param_list; /**< parameter list in string for print */
        EbSvtAv1EncConfiguration
            param_config; /**< configuration of parameter set */
        TestConfig() {
            name = "unnamed";
            description = "";
            use_parser = false;
            is_death = false;
            result = "";
            config_line = "";
            param_list.clear();
            memset(&param_config, 0xFF, sizeof(param_config));
        }
        void print_param() const {
            printf("Name: %s\n", name.c_str());
            printf("Description: %s\n", description.c_str());
            printf("Test Vector: %s\n", test_vector.c_str());
            printf("Use Parser: %s\n", use_parser ? "YES" : "NO");
            printf("Expect Result: %s\n", result.c_str());
            printf("Parameters List:\n");
            for (string param : param_list)
                printf("%s\n", param.c_str());
        }
    } TestConfig;

    typedef void (*set_param)(const string, TestConfig*);

  public:
    ConfigReader(const string& config_path)
        : config_file_(config_path, std::ios::in) {
        config_vec_.clear();
        const TestConfig* config = nullptr;
        while (config = read_config()) {
            config_vec_.push_back(config);
        }
        config_file_.close();
    }
    ~ConfigReader() {
        while (config_vec_.size()) {
            delete config_vec_.back();
            config_vec_.pop_back();
        }
    }
    const int get_config_count() {
        return config_vec_.size();
    }
    const TestConfig* get_config(const int index) {
        if (index >= config_vec_.size())
            return nullptr;
        return config_vec_.at(index);
    }

  private:
    static void set_name(const string value, TestConfig* config) {
        config->name = std::move(value);
    }
    static void set_description(const string value, TestConfig* config) {
        config->description = std::move(value);
    }
    static void set_test_vector(const string value, TestConfig* config) {
        config->test_vector = std::move(value);
    }
    static void set_use_parser(const string value, TestConfig* config) {
        config->use_parser = (std::stoul(value) == 1);
    }
    static void set_is_death(const string value, TestConfig* config) {
        config->is_death = (std::stoul(value) == 1);
    }
    static void set_result(const string value, TestConfig* config) {
        config->result = std::move(value);
    }
    DEFINE_SET_PARAM(enc_mode);
    DEFINE_SET_PARAM(intra_period_length);
    DEFINE_SET_PARAM(intra_refresh_type);
    DEFINE_SET_PARAM(hierarchical_levels);
    DEFINE_SET_PARAM(pred_structure);
    DEFINE_SET_PARAM(base_layer_switch_mode);
    DEFINE_SET_PARAM(source_width);
    DEFINE_SET_PARAM(source_height);
    DEFINE_SET_PARAM(frame_rate);
    DEFINE_SET_PARAM(frame_rate_numerator);
    DEFINE_SET_PARAM(frame_rate_denominator);
    DEFINE_SET_PARAM(encoder_bit_depth);
    static void set_encoder_color_format(const string value,
                                         TestConfig* config) {
        config->param_config.encoder_color_format =
            (EbColorFormat)std::stoul(value);
    }
    DEFINE_SET_PARAM(compressed_ten_bit_format);
    DEFINE_SET_PARAM(frames_to_be_encoded);
    DEFINE_SET_PARAM(improve_sharpness);
    DEFINE_SET_PARAM(sb_sz);
    DEFINE_SET_PARAM(super_block_size);
    DEFINE_SET_PARAM(partition_depth);
    DEFINE_SET_PARAM(qp);
    DEFINE_SET_PARAM(use_qp_file);
    DEFINE_SET_PARAM(enable_qp_scaling_flag);
    DEFINE_SET_PARAM(disable_dlf_flag);
    DEFINE_SET_PARAM(enable_denoise_flag);
    DEFINE_SET_PARAM(film_grain_denoise_strength);
    DEFINE_SET_PARAM(enable_warped_motion);
    DEFINE_SET_PARAM(use_default_me_hme);
    DEFINE_SET_PARAM(enable_hme_flag);
    DEFINE_SET_PARAM(ext_block_flag);
    DEFINE_SET_PARAM(in_loop_me_flag);
    DEFINE_SET_PARAM(search_area_width);
    DEFINE_SET_PARAM(search_area_height);
    DEFINE_SET_PARAM(constrained_intra);
    DEFINE_SET_PARAM(rate_control_mode);
    DEFINE_SET_PARAM(scene_change_detection);
    DEFINE_SET_PARAM(look_ahead_distance);
    DEFINE_SET_PARAM(target_bit_rate);
    DEFINE_SET_PARAM(max_qp_allowed);
    DEFINE_SET_PARAM(min_qp_allowed);
    DEFINE_SET_PARAM(high_dynamic_range_input);
    DEFINE_SET_PARAM(profile);
    DEFINE_SET_PARAM(tier);
    DEFINE_SET_PARAM(level);
    DEFINE_SET_PARAM(asm_type);
    DEFINE_SET_PARAM(channel_id);
    DEFINE_SET_PARAM(active_channel_count);
    DEFINE_SET_PARAM(speed_control_flag);
    DEFINE_SET_PARAM(injector_frame_rate);
    DEFINE_SET_PARAM(logical_processors);
    DEFINE_SET_PARAM(target_socket);
    DEFINE_SET_PARAM(recon_enabled);
#if TILES
    DEFINE_SET_PARAM(tile_columns);
    DEFINE_SET_PARAM(tile_rows);
#endif

  private:
    const TestConfig* read_config() {
        string line;
        if (std::getline(config_file_, line)) {
            /* first trim the string */
            line.erase(std::remove_if(line.begin(), line.end(), std::isspace),
                       line.end());
            TestConfig* new_config = new TestConfig();
            if (new_config) {
                new_config->config_line = std::move(line);
                std::stringstream ss(new_config->config_line);
                string param_str;
                while (std::getline(ss, param_str, ';')) {
                    string::size_type pos = param_str.find("=", 0);
                    if (pos == string::npos) {
                        printf("error: not find parameter pattern!\n");
                        continue;
                    }
                    const string param_name = param_str.substr(0, pos);
                    const string param_value =
                        param_str.substr(pos + 1, param_str.size() - pos);
                    std::map<const string, set_param>::const_iterator it =
                        setting_map_.find(param_name);
                    if (it != setting_map_.end()) {
                        it->second(param_value, new_config);
                    }
                    new_config->param_list.push_back(std::move(param_str));
                }
                if (new_config->param_list.size() == 0) {
                    printf("this is an empty configuration! skip it!\n");
                }
                return new_config;
            }
        }
        return nullptr;
    }

  private:
    static const std::map<const string, set_param> setting_map_;
    std::ifstream config_file_;
    vector<const TestConfig*> config_vec_;
};

const std::map<const string, ConfigReader::set_param>
    ConfigReader::setting_map_ = {
        MAP_GEN(name),
        MAP_GEN(description),
        MAP_GEN(test_vector),
        MAP_GEN(use_parser),
        MAP_GEN(is_death),
        MAP_GEN(result),
        MAP_GEN(enc_mode),
        MAP_GEN(intra_period_length),
        MAP_GEN(intra_refresh_type),
        MAP_GEN(hierarchical_levels),
        MAP_GEN(pred_structure),
        MAP_GEN(base_layer_switch_mode),
        MAP_GEN(source_width),
        MAP_GEN(source_height),
        MAP_GEN(frame_rate),
        MAP_GEN(frame_rate_numerator),
        MAP_GEN(frame_rate_denominator),
        MAP_GEN(encoder_bit_depth),
        MAP_GEN(encoder_color_format),
        MAP_GEN(compressed_ten_bit_format),
        MAP_GEN(frames_to_be_encoded),
        MAP_GEN(improve_sharpness),
        MAP_GEN(sb_sz),
        MAP_GEN(super_block_size),
        MAP_GEN(partition_depth),
        MAP_GEN(qp),
        MAP_GEN(use_qp_file),
        MAP_GEN(enable_qp_scaling_flag),
        MAP_GEN(disable_dlf_flag),
        MAP_GEN(enable_denoise_flag),
        MAP_GEN(film_grain_denoise_strength),
        MAP_GEN(enable_warped_motion),
        MAP_GEN(use_default_me_hme),
        MAP_GEN(enable_hme_flag),
        MAP_GEN(ext_block_flag),
        MAP_GEN(in_loop_me_flag),
        MAP_GEN(search_area_width),
        MAP_GEN(search_area_height),
        MAP_GEN(constrained_intra),
        MAP_GEN(rate_control_mode),
        MAP_GEN(scene_change_detection),
        MAP_GEN(look_ahead_distance),
        MAP_GEN(target_bit_rate),
        MAP_GEN(max_qp_allowed),
        MAP_GEN(min_qp_allowed),
        MAP_GEN(high_dynamic_range_input),
        MAP_GEN(profile),
        MAP_GEN(tier),
        MAP_GEN(level),
        MAP_GEN(asm_type),
        MAP_GEN(channel_id),
        MAP_GEN(active_channel_count),
        MAP_GEN(speed_control_flag),
        MAP_GEN(injector_frame_rate),
        MAP_GEN(logical_processors),
        MAP_GEN(target_socket),
        MAP_GEN(recon_enabled),
#if TILES
        MAP_GEN(tile_columns),
        MAP_GEN(tile_rows),
#endif
};

// #define TEST_MYSELF

#ifdef TEST_MYSELF
#include "gtest/gtest.h"

TEST(ConfigReaderSelfTest, test_myself) {
    ConfigReader* reader = new ConfigReader("./config_test");
    ASSERT_NE(reader, nullptr);
    const ConfigReader::TestConfig* config = nullptr;
    for (int i = 0; i < reader->get_config_count(); i++) {
        config = reader->get_config(i);
        ASSERT_NE(config, nullptr);
        printf("------------- Configuration %d --------------\n", i);
        printf("Configuration Line: %s\n", config->config_line.c_str());
        config->print_param();
        printf("-------------      End %d      --------------\n", i);
    }
}

#endif  // TEST_MYSELF

#endif  // !_CONFIG_READER_H_
