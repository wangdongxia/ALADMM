#include <fstream>
#include <sstream>

#include "logging/simple_logging.h"
#include "data/sparse_dataset.h"
#include "other/string_util.h"

SparseDataset::SparseDataset(const std::string &file_path) : file_path_(file_path), all_feature_num_(0),
                                                             sample_num_(0), dimension_(0) {
    std::ifstream file_reader(file_path_);
    if (file_reader.fail()) {
        LOG(FATAL) << "打开文件失败，文件名：" << file_path_;
    }
    /* 首先完整读取一遍数据集，检查是否有格式错误并且确定特征数量总和，方便创建数组 */
    int label;
    int line_num = 0;
    std::string temp_str;
    std::string line;
    std::stringstream line_reader;
    while (std::getline(file_reader, line)) {
        ++line_num;
        /* 去除首尾空格 */
        Trim(line);
        //line.erase(0, line.find_first_not_of(' '));
        //line.erase(line.find_last_not_of(' ') + 1);
        if (line.empty()) {
            LOG(WARNING) << "出现空行，文件名：" << file_path_ << "，行号：" << line_num;
            continue;
        } else {
            line_reader.clear();
            line_reader.str(line);
            if (!(line_reader >> label) || !(label == 1 || label == -1 || label == 0)) {
                LOG(FATAL) << "读取标签失败，文件名：" << file_path_ << "，行号：" << line_num;
            }
            /* 读取一条样本的所有特征，为键值对格式 */
            while (line_reader >> temp_str) {
                /* 检查格式是否正确 */
                size_t pos = temp_str.find_first_of(':');
                if (pos == std::string::npos || pos == 0 || pos == temp_str.size() - 1) {
                    LOG(FATAL) << "读取特征失败，文件名：" << file_path_ << "，行号：" << line_num;
                }
                std::string index = temp_str.substr(0, pos);
                std::string value = temp_str.substr(pos + 1);
                size_t idx1;
                size_t idx2;
                /* 尝试进行类型转换，检查是否有错误 */
                std::stoi(index, &idx1);
                std::stod(value, &idx2);
                if (idx1 != index.size() || idx2 != value.size()) {
                    LOG(FATAL) << "读取特征失败，文件名：" << file_path_ << "，行号：" << line_num;
                }
                /* 特征索引是从1开始的，这里进行检查 */
                int temp_index = std::stoi(index, &idx1);
                CHECK(temp_index > 0);
                ++all_feature_num_;
            }
            // 每条样本以一个index=-1的Feature表示结束
            ++all_feature_num_;
            ++sample_num_;
        }
    }

    label_list_ = new int[sample_num_];
    sample_list_ = new Feature *[sample_num_];
    sample_space_ = new Feature[all_feature_num_];

    // 让文件读取指针归位
    file_reader.clear();
    file_reader.seekg(0, std::ifstream::beg);
    int i = 0, j = 0;
    while (std::getline(file_reader, line)) {
        Trim(line);
        //line.erase(0, line.find_first_not_of(' '));
        //line.erase(line.find_last_not_of(' ') + 1);
        if (line.empty()) {
            continue;
        } else {
            // 一条样本的起始Feature的地址
            sample_list_[i] = &sample_space_[j];
            line_reader.clear();
            line_reader.str(line);
            line_reader >> label_list_[i];  
            if (label_list_[i] != 1) {           //将标签统一为{1,-1}
                label_list_[i] = -1;
            }
            while (line_reader >> temp_str) {
                size_t pos = temp_str.find_first_of(':');
                std::string index = temp_str.substr(0, pos);
                std::string value = temp_str.substr(pos + 1);
                /* 数据集中的特征索引是从1开始的，因此在这边进行修正为从0开始，方便与数组进行计算 */
                sample_space_[j].index = atoi(index.c_str()) - 1;
                sample_space_[j].value = atof(value.c_str());
                dimension_ = sample_space_[j].index > dimension_ ? sample_space_[j].index : dimension_;
                ++j;
            }
            // 一条样本以index=-1的Feature结尾
            sample_space_[j++].index = -1;
            ++i;
        }
    }
    if (i != sample_num_ || j != all_feature_num_) {
        LOG(FATAL) << "读取数据出错，文件名：" << file_path;
    }
    file_reader.close();
}

SparseDataset::~SparseDataset() {
    delete[] label_list_;
    delete[] sample_list_;
    delete[] sample_space_;
}

int SparseDataset::GetLabel(int n) {
    if (n >= sample_num_) {
        LOG(FATAL) << "读取数据出错，文件名：" << file_path_ << "，总样本数：" << sample_num_ << "，要读取的样本索引：" << n;
    }
    return label_list_[n];
}

const Feature *SparseDataset::GetSample(int n) {
    if (n >= sample_num_) {
        LOG(FATAL) << "读取数据出错，文件名：" << file_path_ << "，总样本数：" << sample_num_ << "，要读取的样本索引：" << n;
    }
    return sample_list_[n];
}
