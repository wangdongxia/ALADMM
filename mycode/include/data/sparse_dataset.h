#ifndef UTILS_SPARSE_DATASET_H
#define UTILS_SPARSE_DATASET_H

#include <string>

struct Feature {
    int index;
    double value;
};

/* 二分类稀疏数据集读取类 */
class SparseDataset {
public:
    explicit SparseDataset(const std::string &file_path);

    ~SparseDataset();

    int GetDimension() { return dimension_; }

    int GetSampleNumber() { return sample_num_; }

    const Feature *GetSample(int n);

    int GetLabel(int n);
private:
    const std::string file_path_;
    int sample_num_;
    /* dimension_是当前数据子集的样本中出现的最大索引，并不代表整个数据集的维度 */
    int dimension_;
    /* 为了高效读取，整个稀疏数据集的特征被存储在一个数组中，这个值表示所有样本特征数量的总和 */
    int all_feature_num_;
    int *label_list_;
    /* 为了方便使用，储存每条样本的起始地址，无需储存结束地址，index=-1代表样本结尾 */
    Feature **sample_list_;
    /* 用来存储所有特征 */
    Feature *sample_space_;
};

#endif //SPARSE_DATASET_H
