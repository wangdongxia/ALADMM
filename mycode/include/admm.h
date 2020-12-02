#ifndef RADADMM_ADMM_H
#define RADADMM_ADMM_H

#include <string>
#include <mutex>

#include "sparse_allreduce.h"
#include "data/sparse_dataset.h"
#include "optimizer/lr_tron_optimizer.h"

class ADMM {
public:
    ADMM(int dimension, int min_barrier, int max_delay, int max_iterations, int interval, double rho, double l2reg, double l1reg,
         double ABSTOL, double RELTOL, std::string train_data_path, std::string test_data_path,int type,double theta);

    ~ADMM();

    int GetID() { return id_; }

    void Run();

private:
    double theta_;
	int filter_type_;
    bool stopped_;
    bool triggered_;
    int k_;
    int id_;
    int worker_number_;
    int dimension_;
    int max_iterations_;
    double rho_;
    double l2reg_;
	double l1reg_;
    double calculate_time_;
    double ABSTOL_, RELTOL_;
    double *x_, *y_, *z_, *w_;
	double *delta_w_;
	double *subx_z_;
    double *old_x_, *old_y_, *old_z_;
    double *temp_z_;
    timeval start_time_, end_time_;
    std::mutex m_;
    SparseDataset *train_data_;
    SparseDataset *test_data_;
    LRTronOptimizer *optimizer_;
    spar::Coordinator *coordinator_;
    spar::Communicator<spar::SumOperator, double> *communicator_;

    static void PostAsynAllreduce(double *buffer, int dimension, void *parameter);

    static void SupplyAllreduceData(double *buffer, int dimension, void *parameter);
};

#endif //RADADMM_ADMM_H
