#ifndef UTILS_TRONOPTIMIZER_H
#define UTILS_TRONOPTIMIZER_H

#include "optimizer/optimizer.h"
#include "data/sparse_dataset.h"

class LRTronOptimizer : public Optimizer {
public:
    LRTronOptimizer(const double *y, const double *z, int dimension, double rho, int max_iterations,
                    double epsilon, double cg_epsilon, SparseDataset *dataset);

    ~LRTronOptimizer();

    void Optimize(double *x);

    void SetDimension(int dimension) { dimension_ = dimension; }

    void SetMaxIterations(int max_iterations) { max_iterations_ = max_iterations; }

    void SetCgEpsilon(double cg_epsilon) { cg_epsilon_ = cg_epsilon; }

    void SetEpsilon(double epsilon) { epsilon_ = epsilon; }

    void SetY(const double *y) { y_ = y; }

    void SetZ(const double *z) { z_ = z; }

    void SetDataset(SparseDataset *dataset);

private:
    int dimension_;
    int max_iterations_;
    double epsilon_;
    double tron_epsilon_;
    double cg_epsilon_;
    double rho_;
    double *D_;
    const double *y_;
    const double *z_;
    SparseDataset *dataset_;

    int trpcg(double delta, double *g, double *M, double *s, double *r, bool *reach_boundary);

    double function_value(const double *x);

    void gradient(const double *x, double *g);

    void get_diag_preconditioner(double *M);

    void Hv(double *s, double *Hs);
};

#endif
