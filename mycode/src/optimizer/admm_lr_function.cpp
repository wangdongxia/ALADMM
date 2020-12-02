#include "optimizer/admm_lr_function.h"
#include "math/simple_algebra.h"

double AdmmLRFunction::Evaluate(const double *x) {
    double sum = 0.0;
    for (int i = 0; i < dimension_; ++i) {
        double temp = x[i] - z_[i] + y_[i] / rho_;
        sum += temp * temp;
    }
    sum *= (rho_ / 2.0);
    int sample_num = dataset_->GetSampleNumber();
    for (int i = 0; i < sample_num; ++i) {
        sum += std::log(1 + exp(-dataset_->GetLabel(i) * Dot(x, dataset_->GetSample(i))));
    }
    return sum;
}

void AdmmLRFunction::Gradient(const double *x, double *g) {
    for (int i = 0; i < dimension_; ++i) {
        g[i] = y_[i] + rho_ * (x[i] - z_[i]);
    }
    int sample_num = dataset_->GetSampleNumber();
    for (int i = 0; i < sample_num; ++i) {
        double d = (Sigmoid(dataset_->GetLabel(i) * Dot(x, dataset_->GetSample(i))) - 1) * dataset_->GetLabel(i);
        const Feature *sample = dataset_->GetSample(i);
        while (sample->index != -1) {
            g[sample->index] += (sample->value * d);
            ++sample;
        }
    }
}
