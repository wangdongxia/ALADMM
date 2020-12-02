#ifndef UTILS_ADMM_LR_FUNCTION_H
#define UTILS_ADMM_LR_FUNCTION_H

#include "optimizer/differentiable_function.h"
#include "data/sparse_dataset.h"

/* ADMM算法下LR模型的损失函数 */
class AdmmLRFunction : public DifferentiableFunction {
public:
    AdmmLRFunction(const double *y, const double *z, int dimension, double rho,
                   SparseDataset *dataset) : y_(y), z_(z), rho_(rho), dimension_(dimension), dataset_(dataset) {}

    double Evaluate(const double *x);

    void Gradient(const double *x, double *g);

    void SetRho(double rho) { rho_ = rho; }

    void SetDataset(SparseDataset *dataset) { dataset_ = dataset; }

    void SetY(const double *y) { y_ = y; }

    void SetZ(const double *z) { z_ = z; }

    void SetDimension(int dimension) { dimension_ = dimension; }

private:
    int dimension_;
    double rho_;
    const double *y_;
    const double *z_;
    SparseDataset *dataset_;
};

#endif
