#ifndef SPARSEALLREDUCE_COMMUNICATOR_H
#define SPARSEALLREDUCE_COMMUNICATOR_H

#include <thread>
#include <vector>
#include <condition_variable>
#include <mutex>
#include <map>
#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>

#include "logging/simple_logging.h"
#include "internal/common.h"
#include "internal/group_allreduce.h"
#include "internal/ring_allreduce.h"
#include "internal/split_sparse_allreduce.h"
#include "internal/simple_broadcast.h"
#include "internal/communicator_singleton.h"

namespace spar {

template<class O, class T>
class Communicator {
public:
    static Communicator *GetInstance();

    void Run();

    int AsynAllreduce(T *buffer);

    template<class O2, class T2>
    void SyncAllreduce(T2 *buffer, int count);

    template<class T2>
    void SyncBroadcast(T2 *buffer, int count, int root);

    void WaitUntilCommunicatorStop();

    int GetID() { return id_; }

    int GetWorkerNumber() { return worker_number_; }

    void SetCallbackParameter(void *parameter) { callback_parameter_ = parameter; }

    void SetDimension(int dimension) { dimension_ = dimension; }
	void SetFilterType(int type) { filter_type_ = type; }

    void SetPreAsynAllreduce(void PreAsynAllreduce(T *, int, void *)) { PreAsynAllreduce_ = PreAsynAllreduce; }

    void SetPostAsynAllreduce(void PostAsynAllreduce(T *, int, void *)) { PostAsynAllreduce_ = PostAsynAllreduce; }

    void SetSupplyAllreduceData(void SupplyAllreduceData(T *, int, void *)) {
        SupplyAllreduceData_ = SupplyAllreduceData;
    }

    void SetPreTerminate(void PreTerminate(void *)) { PreTerminate_ = PreTerminate; }

private:
    int filter_type_;
    int id_;
    int coordinator_id_;
    int worker_number_;
    int dimension_;
    double comm_time_;
    
    T *user_buffer_;
     
    void *callback_parameter_;
     
    std::vector<int> member_list_;
     
    std::vector<int> leader_list_;
     
    std::vector<int> worker_list_;
    CalculatorState calculator_state_;
    CommunicatorState communicator_state_;
    std::mutex m_;
    std::condition_variable cv_;

    void (*PreAsynAllreduce_)(T *, int, void *);

    void (*PostAsynAllreduce_)(T *, int, void *);

    void (*SupplyAllreduceData_)(T *, int, void *);

    void (*PreTerminate_)(void *);

    Communicator();

    void CreateGroup();

    void ThreadFunction();

    void Terminate();

    bool ContainSelf(int *ready_worker_list, int ready_worker_num);
};

template<class O, class T>
Communicator<O, T> *Communicator<O, T>::GetInstance() {
    
    static Communicator communicator;
    CommunicatorSingleton::Set(&communicator);
    return &communicator;
}

template<class O, class T>
Communicator<O, T>::Communicator() :dimension_(-1), PreAsynAllreduce_(NULL), PostAsynAllreduce_(NULL), comm_time_(0),
                                    SupplyAllreduceData_(NULL), PreTerminate_(NULL), callback_parameter_(NULL) {
    int flag;
    MPI_Initialized(&flag);
    CHECK_EQ(flag, 1) << "pleas start the Communicator thread after spar::Init()";
    MPI_Comm_rank(MPI_COMM_WORLD, &id_);
    MPI_Comm_size(MPI_COMM_WORLD, &worker_number_);
     
    --worker_number_;
    
    coordinator_id_ = worker_number_;
    calculator_state_ = CalculatorState::kCalculating;
    communicator_state_ = CommunicatorState::kWaitToRun;
}

template<class O, class T>
void Communicator<O, T>::WaitUntilCommunicatorStop() {
    m_.lock();
    if (calculator_state_ != CalculatorState::kStopped) {
        calculator_state_ = CalculatorState::kWantToStop;
    }
    m_.unlock();
    cv_.notify_all();
    while (communicator_state_ != CommunicatorState::kStopped) {
        sleep(1);
    }
}

template<class O, class T>
template<class O2, class T2>
void Communicator<O, T>::SyncAllreduce(T2 *buffer, int count) {
	 
	if(filter_type_==3)
	{
		SplitAlgather<O2>(buffer, count, id_, worker_list_, MPI_COMM_WORLD);
		 
	}else if(filter_type_==2){
		HierSplitAlgather<O2>(buffer, count, id_, worker_list_, MPI_COMM_WORLD);
	}else if(filter_type_==1){
		GroupRingAllreduce<O2>(buffer, count, id_, member_list_, leader_list_, MPI_COMM_WORLD);
	}
	else{ 
		 RingAllreduce<O2>(buffer, count, id_, worker_list_, MPI_COMM_WORLD);
    }
}

template<class O, class T>
template<class T2>
void Communicator<O, T>::SyncBroadcast(T2 *buffer, int count, int root) {
    SimpleBroadcast(buffer, count, root, id_, worker_list_, MPI_COMM_WORLD);
}

template<class O, class T>
int Communicator<O, T>::AsynAllreduce(T *buffer) {
    std::unique_lock<std::mutex> lck(m_);
    
    if (calculator_state_ == CalculatorState::kWantToStop ||
        calculator_state_ == CalculatorState::kStopped) {
        return 1;
    }
    
    user_buffer_ = buffer;
     
    calculator_state_ = CalculatorState::kCalculationComplete;
    while (calculator_state_ != CalculatorState::kAllreduceComplete &&
           calculator_state_ != CalculatorState::kWantToStop &&
           calculator_state_ != CalculatorState::kStopped) {
        cv_.wait(lck);
    }
     
    if (calculator_state_ == CalculatorState::kWantToStop ||
        calculator_state_ == CalculatorState::kStopped) {
        return 1;
    } else {
        calculator_state_ = CalculatorState::kCalculating;
        return 0;
    }
}

template<class O, class T>
void Communicator<O, T>::Run() {
    CHECK(dimension_ != -1 && SupplyAllreduceData_ != NULL) << "please set parameters";
    std::unique_lock<std::mutex> lck(m_);
    
    if (communicator_state_ == CommunicatorState::kWaitToRun) {
        
        std::thread t(&Communicator::ThreadFunction, this);
        t.detach();
        communicator_state_ = CommunicatorState::kRunning;
    }
}

template<class O, class T>
void Communicator<O, T>::ThreadFunction() {
    
	 
    CreateGroup();
    int flag;
    bool stopped = false;
    MPI_Status status;
    int ready_worker_num;
    int *ready_worker_list = new int[worker_number_];
    T *temp_buffer = new T[dimension_];
    timeval start_time, end_time;
    while (!stopped) {
        m_.lock();
         
        if (calculator_state_ == CalculatorState::kCalculationComplete) {
            MPI_Request temp;
            MPI_Isend(NULL, 0, MPI_INT, coordinator_id_, MessageType::kAllreduceApplication, MPI_COMM_WORLD,
                      &temp);
            MPI_Request_free(&temp);
            calculator_state_ = CalculatorState::kWaitingAllreduceResult;
        }
        
        if (calculator_state_ == CalculatorState::kWantToStop) {
            MPI_Request temp;
            MPI_Isend(NULL, 0, MPI_INT, coordinator_id_, MessageType::kTerminateApplication, MPI_COMM_WORLD,
                      &temp);
            MPI_Request_free(&temp);
            calculator_state_ = CalculatorState::kStopped;
        }
        m_.unlock();
        MPI_Iprobe(coordinator_id_, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            switch (status.MPI_TAG) {
                 
                case MessageType::kAllreduceCommand:
                    MPI_Get_count(&status, MPI_INT, &ready_worker_num);
                    
                    MPI_Recv(ready_worker_list, ready_worker_num, MPI_INT, coordinator_id_,
                             MessageType::kAllreduceCommand, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    T *allreduce_buffer;
                    if (ContainSelf(ready_worker_list, ready_worker_num)) {
                        allreduce_buffer = user_buffer_;
                    } else {
                        allreduce_buffer = temp_buffer;
                         
                        SupplyAllreduceData_(allreduce_buffer, dimension_, callback_parameter_);
                    }
                     
                    if (PreAsynAllreduce_) {
                        PreAsynAllreduce_(allreduce_buffer, dimension_, callback_parameter_);
                    }
                    gettimeofday(&start_time, NULL);
                    SyncAllreduce<O, T>(allreduce_buffer, dimension_);
                    gettimeofday(&end_time, NULL);
                    comm_time_ += ((end_time.tv_sec - start_time.tv_sec) +
                                   (end_time.tv_usec - start_time.tv_usec) / 1000000.0);
 
                    if (PostAsynAllreduce_) {
                        PostAsynAllreduce_(allreduce_buffer, dimension_, callback_parameter_);
                    }
                    MPI_Request temp;
                    MPI_Isend(NULL, 0, MPI_INT, coordinator_id_, MessageType::kAllreduceComplete, MPI_COMM_WORLD,
                              &temp);
                    MPI_Request_free(&temp);
                     
                    if (ContainSelf(ready_worker_list, ready_worker_num)) {
                        std::unique_lock<std::mutex> lck(m_);
                        calculator_state_ = CalculatorState::kAllreduceComplete;
                        cv_.notify_all();
                    }
                    break;
                case MessageType::kTerminateCommand:
                    MPI_Recv(NULL, 0, MPI_INT, coordinator_id_, MessageType::kTerminateCommand, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    stopped = true;
                    break;
                default:
                    LOG(FATAL) << "unknown message";
            }
        } else {
            usleep(100000);
        }
    }
    delete[] temp_buffer;
    delete[] ready_worker_list;
    Terminate();
}

template<class O, class T>
void Communicator<O, T>::CreateGroup() {
    char name[50];
    int name_length;
     
    MPI_Get_processor_name(name, &name_length);
    MPI_Send(name, name_length, MPI_CHAR, coordinator_id_, MessageType::kCreateGroup1, MPI_COMM_WORLD);
    MPI_Status status;
    int list_length;
    MPI_Probe(coordinator_id_, MessageType::kCreateGroup1, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_INT, &list_length);
    int *recv_buf = new int[list_length];
    MPI_Recv(recv_buf, list_length, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
    member_list_.assign(recv_buf, recv_buf + list_length);
    delete[] recv_buf;
    
    int leader_id_ = member_list_.front();
    if (id_ == leader_id_) {
        
        MPI_Probe(coordinator_id_, MessageType::kCreateGroup2, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &list_length);
        recv_buf = new int[list_length];
        MPI_Recv(recv_buf, list_length, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        leader_list_.assign(recv_buf, recv_buf + list_length);
        delete[] recv_buf;
    }
    for (int i = 0; i < worker_number_; ++i) {
        worker_list_.push_back(i);
    }
}

template<class O, class T>
void Communicator<O, T>::Terminate() {
    SyncAllreduce<SumOperator>(&comm_time_, 1);
    if (id_ == 0) {
        std::cout << comm_time_ / worker_number_ << std::endl;
    }
    if (PreTerminate_) {
        PreTerminate_(callback_parameter_);
    }
    MPI_Send(NULL, 0, MPI_INT, coordinator_id_, MessageType::kTerminateComplete, MPI_COMM_WORLD);
    m_.lock();
    communicator_state_ = CommunicatorState::kStopped;
    calculator_state_ = CalculatorState::kStopped;
    m_.unlock();
    cv_.notify_all();
    sleep(3);
 
}

template<class O, class T>
bool Communicator<O, T>::ContainSelf(int *ready_worker_list, int ready_worker_num) {
    for (int i = 0; i < ready_worker_num; ++i) {
        if (ready_worker_list[i] == id_) {
            return true;
        } else if (ready_worker_list[i] > id_) {
            return false;
        }
    }
    return false;
}

}

#endif //SPARSEALLREDUCE_COMMUNICATOR_H
