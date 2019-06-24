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
#include "E2eTestVectors.h"

using std::string;
using std::vector;
using namespace svt_av1_e2e_test_vector;

#define MAX_CONFIG_STR_LEN (1024)

#define CFG_LINE_KEY_STR "Configure"

#define MAP_GEN(param_name) \
    { #param_name, ConfigReader::set_##param_name }

class ConfigReader {
  public:
    typedef struct TestConfig : public EncTestSetting {
        string description;     /**< description of test case in string */
        string test_vector_cfg; /**< config file of test vector for this case */
        string config_line;     /**< configuration in a line */
        vector<string> arg_vec; /**< argument list from config_line */

        TestConfig() {
            description = "";
            test_vector_cfg = "";
            config_line = "";
            arg_vec.clear();
        }

        /** generate TestConfig from a EncTestSetting, only for parsing of
         * parameter line */
        TestConfig(const EncTestSetting& src) {
            name = src.name;
            setting = src.setting;
            test_vectors = src.test_vectors;
            description = src.name;
            test_vector_cfg = "";
            config_line = setting.find(CFG_LINE_KEY_STR)->second;
            arg_vec.clear();
            parse_cfg_line();
        }

        void print_param() const {
            printf("Name: %s\n", name.c_str());
            printf("Description: %s\n", description.c_str());
            printf("Test Vector: %s\n", test_vector_cfg.c_str());
            printf("Configuration Line: %s\n", config_line.c_str());
        }

        void parse_cfg_line() {
            /** parse arguments to list */
            std::stringstream ss_config(config_line);
            string arg_str;
            while (std::getline(ss_config, arg_str, ' ')) {
                if (arg_str.length() > 0)
                    arg_vec.push_back(arg_str);
            }
        }

        int get_arg_arr(const char* str_arr[], const int max_count) const {
            int arg_num = static_cast<int>(arg_vec.size());
            if (arg_num > 0) {
                for (int i = 0; i < std::min(arg_num, max_count); i++) {
                    str_arr[i] = arg_vec.at(i).c_str();
                }
            }
            return arg_num;
        }
    } TestConfig;

    typedef void (*set_param)(const string, TestConfig*);

  public:
    ConfigReader(const string& config_path)
        : config_file_(config_path, std::ios::in) {
        config_vec_.clear();
        const TestConfig* config = nullptr;
        while ((config = read_config()) != nullptr) {
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

    int get_config_count() {
        return static_cast<int>(config_vec_.size());
    }

    const TestConfig* get_config(const int index) {
        if (index >= static_cast<int>(config_vec_.size()))
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
        config->test_vector_cfg = std::move(value);
    }

    static void set_config_line(const string value, TestConfig* config) {
        config->config_line = std::move(value);
    }

  private:
    const TestConfig* read_config() {
        do {
            string line;
            if (std::getline(config_file_, line)) {
                /** skip the comments and empty line */
                if (line.size() == 0)
                    continue;
                const char start_char = line.at(0);
                if (start_char == '#')
                    continue;
                TestConfig* new_config = new TestConfig();
                if (new_config) {
                    std::stringstream ss(line);
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
                    }
                    if (new_config->config_line.length() > 0) {
                        new_config->setting = {
                            {CFG_LINE_KEY_STR,
                             new_config->config_line.c_str()}};
                        new_config->parse_cfg_line();
                    }
                    if (new_config->test_vector_cfg.length() > 0) {
                        new_config->test_vectors = generate_vector_from_config(
                            new_config->test_vector_cfg.c_str());
                    }
                    return new_config;
                }
            } else
                break;
        } while (1);

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
        MAP_GEN(config_line),
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
        config->print_param();
        printf("-------------      End %d      --------------\n", i);
    }
}

TEST(ConfigReaderSelfTest, copy_setting) {
    const EncTestSetting setting = {
        "copy_setting_test",
        {{CFG_LINE_KEY_STR, "-i c:/test.yuv -q 30 -o output.ivf"}},
        default_test_vectors};
    const ConfigReader::TestConfig config(setting);
    printf("------------- Configuration  --------------\n");
    config.print_param();
    int max_args = 100;
    const char* args[100] = {0};
    int args_num = config.get_arg_arr(args, max_args);
    ASSERT_GT(args_num, 0);
    ASSERT_LE(args_num, max_args);
    for (size_t i = 0; i < args_num; i++) {
        printf("%s\n", args[i]);
    }
    printf("-------------      End       --------------\n");
}

#endif  // TEST_MYSELF

#endif  // !_CONFIG_READER_H_
