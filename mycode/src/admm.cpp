#include <sys/time.h>

#include "admm.h"
#include "math/simple_algebra.h"

double Predict(const double *z, SparseDataset *test_dataset) {
    int counter = 0;
    int sample_num = test_dataset->GetSampleNumber();
    for (int i = 0; i < sample_num; ++i) {
        double temp = 1.0 / (1 + exp(-1 * Dot(z, test_dataset->GetSample(i))));
        if (test_dataset->GetLabel(i) == 1 && temp >= 0.5) {
            ++counter;
        }
        if (test_dataset->GetLabel(i) == -1 && temp < 0.5) {
            ++counter;
        }
    }
    return counter * 100.0 / sample_num;
}

double ObjectiveValue(const double *z, SparseDataset *test_dataset) {
    double sum = 0;
    int sample_num = test_dataset->GetSampleNumber();
    for (int i = 0; i < sample_num; ++i) {
        sum += std::log(1 + std::exp(-test_dataset->GetLabel(i) * Dot(z, test_dataset->GetSample(i))));
    }
    return sum/sample_num;
}

ADMM::ADMM(int dimension, int min_barrier, int max_delay, int max_iterations, int interval, double rho,
           double l2reg, double l1reg, double ABSTOL, double RELTOL, std::string train_data_path,
           std::string test_data_path,int type,double theta) : dimension_(dimension), max_iterations_(max_iterations), rho_(rho),
                                         l2reg_(l2reg), l1reg_(l1reg), ABSTOL_(ABSTOL), RELTOL_(RELTOL),filter_type_(type),theta_(theta),
                                         triggered_(false), stopped_(false), k_(0), id_(-1), worker_number_(-1),
                                         calculate_time_(0), optimizer_(NULL), train_data_(NULL), test_data_(NULL),
                                         x_(NULL), y_(NULL), z_(NULL), w_(NULL),delta_w_(NULL),subx_z_(NULL),
                                         old_x_(NULL), old_y_(NULL), old_z_(NULL), temp_z_(NULL) {
    spar::Init(NULL, NULL);
    if (!spar::IsWorker()) {
        coordinator_ = spar::Coordinator::GetInstance();
        coordinator_->SetMaxDelay(max_delay);
        coordinator_->SetMinBarrier(min_barrier);
        coordinator_->SetInterval(interval);
    } else {
        communicator_ = spar::Communicator<spar::SumOperator, double>::GetInstance();
        id_ = communicator_->GetID();
        worker_number_ = communicator_->GetWorkerNumber();
        communicator_->SetCallbackParameter(this);
        communicator_->SetDimension(dimension_);
		communicator_->SetFilterType(filter_type_);
        communicator_->SetPostAsynAllreduce(&ADMM::PostAsynAllreduce);
        communicator_->SetSupplyAllreduceData(&ADMM::SupplyAllreduceData);
        x_ = new double[dimension_];
        y_ = new double[dimension_];
        z_ = new double[dimension_];
        w_ = new double[dimension_];
		delta_w_ = new double[dimension_];
		subx_z_ = new double[dimension_];
        old_x_ = new double[dimension_];
        old_y_ = new double[dimension_];
        old_z_ = new double[dimension_];
        temp_z_ = new double[dimension_];
        char real_path[50];
        sprintf(real_path, train_data_path.c_str(), id_);
        train_data_ = new SparseDataset(real_path);
        if (id_ == 0) {
            test_data_ = new SparseDataset(test_data_path);
        }
        optimizer_ = new LRTronOptimizer(y_, z_, dimension_, rho_, 1, 1e-4, 0.1, train_data_);
    }
}

ADMM::~ADMM() {
    delete[] x_;
    delete[] y_;
    delete[] z_;
    delete[] w_;
	delete[] delta_w_;
	delete[] subx_z_;
    delete[] old_x_;
    delete[] old_y_;
    delete[] old_z_;
    delete[] temp_z_;
    delete optimizer_;
    delete train_data_;
    delete test_data_;
}

void ADMM::Run() {
    if (id_ == -1) {
        coordinator_->Run();
    } else {
        FillZero(x_, dimension_);
        FillZero(y_, dimension_);
        FillZero(z_, dimension_);
        FillZero(w_, dimension_);
		FillZero(delta_w_, dimension_);
		FillZero(subx_z_, dimension_);
        FillZero(old_x_, dimension_);
        FillZero(old_y_, dimension_);
        FillZero(old_z_, dimension_);
        FillZero(temp_z_, dimension_);
        communicator_->Run();
        gettimeofday(&start_time_, NULL);
        timeval cal_start_time, cal_end_time;
        if (id_ == 0) {
            printf("%3s %10s %10s %10s %10s %10s %10s %10s\n", "#", "r norm", "esp_pri", "s norm", "esp_dual",
                   "obj_val", "accuracy", "time");
        }
		int count=0;
        while (true) {
            gettimeofday(&cal_start_time, NULL);
            optimizer_->Optimize(x_);
           /*  for (int i = 0; i < dimension_; ++i) {
                w_[i] = rho_ * x_[i] + y_[i];
            }  */
			 
			if((count>0)&&filter_type_>1)
			{
				 
				for (int i = 0; i < dimension_; ++i) {
					   subx_z_[i]+=rho_ * (x_[i] - old_x_[i]) + (y_[i] - old_y_[i]);
                     if (std::fabs(subx_z_[i]) > theta_*(1.0/sqrt(count+1))) 
					 {   //过滤参数
                        delta_w_[i]=subx_z_[i];
					    subx_z_[i]=0;
						 
                    } else {
                        delta_w_[i] = 0;
						 
                    }   
					/*if (std::abs(y_[i] - old_y_[i]) < theta_ * rho_) {
                        delta_w_[i] = 0;
                    } else {
                        delta_w_[i] = rho_ * (x_[i] - old_x_[i]) + (y_[i] - old_y_[i]);
                    }*/
				}
				 
			}else{
				for (int i = 0; i < dimension_; ++i)
				{
					delta_w_[i]= rho_ * (x_[i]-old_x_[i]) + y_[i]-old_y_[i];
				}
			}
			count++;
            gettimeofday(&cal_end_time, NULL);
            m_.lock();
            calculate_time_ += ((cal_end_time.tv_sec - cal_start_time.tv_sec) +
                                (cal_end_time.tv_usec - cal_start_time.tv_usec) / 1000000.0);
            m_.unlock();
            if (communicator_->AsynAllreduce(delta_w_) == 0) {
				if (stopped_) break;
                for (int i = 0; i < dimension_; ++i) {
                    y_[i] += rho_ * (x_[i] - z_[i]);
                }
            } else {
                break;
            }
        }
        communicator_->WaitUntilCommunicatorStop();
    }
    spar::Finalize();
}

void ADMM::SupplyAllreduceData(double *buffer, int dimension, void *parameter) {
    ADMM *admm = static_cast<ADMM *>(parameter);
    admm->triggered_ = true;
    for (int i = 0; i < dimension; ++i) {
      //  buffer[i] = admm->rho_ * admm->old_x_[i] + admm->old_y_[i];
		  buffer[i]= 0;
    }
}

void ADMM::PostAsynAllreduce(double *buffer, int dimension, void *parameter) {
    ADMM *admm = static_cast<ADMM *>(parameter);
    ++(admm->k_);
	 for (int i = 0; i < dimension; ++i) {
        admm->w_[i] += buffer[i];
    } 
    double *z;
    if (!admm->triggered_) {
        Assign(admm->old_x_, admm->x_, dimension);
        Assign(admm->old_y_, admm->y_, dimension);
        z = admm->z_;
    } else {
        z = admm->temp_z_;
    }
    int worker_number = admm->worker_number_;
    double rho = admm->rho_;
    double l2reg = admm->l2reg_;
    double ABSTOL = admm->ABSTOL_;
    double RELTOL = admm->RELTOL_;
	double l1reg=admm->l1reg_;
    double temp = 1.0/(2 * l2reg + worker_number * rho);
	double t= l1reg*temp;
    for (int i = 0; i < dimension; ++i) {
        z[i] = admm->w_[i] *temp;
		if(l1reg>0)
		{
			if(z[i]>t)
			{
				z[i]-=t;
			}else if(z[i]<-t)
			{
				z[i]+=t;
			}else{
				z[i]=0.0;
			}
		}
    }
    double temp_buffer[3] = {0, 0, 0};
    for (int i = 0; i < dimension; ++i) {
        temp_buffer[0] += admm->old_x_[i] * admm->old_x_[i];
        temp_buffer[1] += admm->old_y_[i] * admm->old_y_[i];
        double temp = admm->old_x_[i] - z[i];
        temp_buffer[2] += temp * temp;
    }
    admm->communicator_->SyncAllreduce<spar::SumOperator>(temp_buffer, 3);
    double nxstack = sqrt(temp_buffer[0]); /* sqrt(sum ||x_i||_2^2) */
    double nystack = sqrt(temp_buffer[1]); /* sqrt(sum ||y_i||_2^2) */
    double prires = sqrt(temp_buffer[2]); /* sqrt(sum ||r_i||_2^2) */
    double z_diff = 0; /* 存放||z_new - z_old||_2^2 */
    double z_norm = 0; /* 存放||z_new||_2^2 */
    for (int i = 0; i < dimension; ++i) {
        double temp = admm->old_z_[i] - z[i];
        z_diff += temp * temp;
        z_norm += z[i] * z[i];
    }
    double dualres = rho * sqrt(worker_number * z_diff);
    double eps_pri = sqrt(dimension * worker_number) * ABSTOL + RELTOL * fmax(nxstack, sqrt(worker_number * z_norm));
    double eps_dual = sqrt(dimension * worker_number) * ABSTOL + RELTOL * nystack;
    if (admm->id_ == 0) {
        gettimeofday(&admm->end_time_, NULL);
        double wait_time = (admm->end_time_.tv_sec - admm->start_time_.tv_sec) +
                           (admm->end_time_.tv_usec - admm->start_time_.tv_usec) / 1000000.0;
        double temp_ov = ObjectiveValue(z, admm->test_data_);
        double temp_ac = Predict(z, admm->test_data_);
        printf("%3d %10.4f %10.4f %10.4f %10.4f %10.6f %10.4f %10.4f,%10.4f\n", admm->k_, prires, eps_pri, dualres,
               eps_dual, temp_ov, temp_ac, wait_time,admm->calculate_time_);
    }
    if (admm->k_ >= admm->max_iterations_) {
//    if (admm->k_ >= admm->max_iterations_ || (prires <= eps_pri && dualres <= eps_dual)) {
        admm->stopped_ = true;
        gettimeofday(&admm->end_time_, NULL);
        double total_time = (admm->end_time_.tv_sec - admm->start_time_.tv_sec) +
                            (admm->end_time_.tv_usec - admm->start_time_.tv_usec) / 1000000.0;
        double temp[] = {total_time, 0};
        admm->m_.lock();
        temp[1] = admm->calculate_time_;
        admm->m_.unlock();
        admm->communicator_->SyncAllreduce<spar::SumOperator>(temp, 2);
        if (admm->id_ == 0) {
            std::cout << "total time:" << temp[0] / admm->worker_number_ << std::endl;
            std::cout << "compute time:" << temp[1] / admm->worker_number_ << std::endl;
            std::cout << "wait time:" << (temp[0] - temp[1]) / admm->worker_number_ << std::endl;
        }
    }
    Assign(admm->old_z_, z, dimension);
    admm->triggered_ = false;
}

