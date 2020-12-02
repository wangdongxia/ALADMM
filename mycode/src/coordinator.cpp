#include <vector>
#include <algorithm>
#include <mpi.h>
#include <unistd.h>

#include "logging/simple_logging.h"
#include "internal/common.h"
#include "coordinator.h"

namespace spar {

Coordinator *Coordinator::GetInstance() {
    //设计成单例模式
    static Coordinator coordinator;
    return &coordinator;
}

Coordinator::Coordinator() : min_barrier_(-1), max_delay_(-1), interval_(-1), k_(0) {
    int flag;
    MPI_Initialized(&flag);
    CHECK_EQ(flag, 1) << "请在调用spar::Init()之后启动Coordinator";
    MPI_Comm_rank(MPI_COMM_WORLD, &id_);
    MPI_Comm_size(MPI_COMM_WORLD, &worker_number_);
    // Worker的数量为总进程数减一
    --worker_number_;
    // Coordinator的id应该为总进程数减一
    CHECK_EQ(id_, worker_number_) << "不要在Worker节点上启动Coordinator";
    for (int i = 0; i < worker_number_; ++i) {
        worker_delay_[i] = 0;
    }
}

void Coordinator::Run() {
    CHECK(min_barrier_ != -1 && max_delay_ != -1) << "请设置必要的参数！";
    //首先进行进程分组
    CreateGroup();
    bool stopped = false;
    MPI_Status status;
    while (!stopped) {
        MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch (status.MPI_TAG) {
            case MessageType::kAllreduceApplication:
                ready_worker_list_.push_back(status.MPI_SOURCE);
                worker_delay_[status.MPI_SOURCE] = -1;
                //判断两个条件是否满足，满足的话命令所有Worker进行Allreduce
                if (IsSatisfied()) {
                    std::sort(ready_worker_list_.begin(), ready_worker_list_.end());
                    MPI_Request temp_request;
                    for (int i = 0; i < worker_number_; ++i) {
                        MPI_Isend(ready_worker_list_.data(), ready_worker_list_.size(), MPI_INT, i,
                                  MessageType::kAllreduceCommand, MPI_COMM_WORLD, &temp_request);
                        MPI_Request_free(&temp_request);
                    }
                    //等待所有Worker完成Allreduce
                    std::vector<MPI_Request> requests;
                    for (int i = 0; i < worker_number_; ++i) {
                        requests.emplace_back();
                        MPI_Irecv(NULL, 0, MPI_INT, i, MessageType::kAllreduceComplete, MPI_COMM_WORLD,
                                  &requests.back());
                    }
                    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
                    for (auto it = worker_delay_.begin(); it != worker_delay_.end(); ++it) {
                        ++it->second;
                    }
                    ready_worker_list_.clear();
                    ++k_;
                }
                break;
            case MessageType::kTerminateApplication:
                stopped = true;
                break;
            default:
                LOG(FATAL) << "未知的消息类型";
        }
    }
    Terminate();
}

void Coordinator::CreateGroup() {
    int count = 0;
    char name[50];
    int name_length;
    std::string hostname;
    std::vector<int> leader_list;
    // 主机名->{id1, id2, id3,...}，其中第一个进程为leader
    std::map<std::string, std::vector<int>> groups;
    MPI_Status status;
    // 等待所有进程将自己所在的主机名发送过来以便进行分组
    while (count < worker_number_) {
        MPI_Probe(MPI_ANY_SOURCE, MessageType::kCreateGroup1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &name_length);
        MPI_Recv(name, name_length, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        hostname.assign(name, name_length);
        groups[hostname].push_back(status.MPI_SOURCE);
        ++count;
    }
    for (auto map_it = groups.begin(); map_it != groups.end(); ++map_it) {
        // 将每个组进行排序，id最小的作为leader
        std::sort(map_it->second.begin(), map_it->second.end());
        leader_list.push_back(map_it->second.front());
        for (auto list_it = map_it->second.begin(); list_it != map_it->second.end(); ++list_it) {
            MPI_Send(map_it->second.data(), map_it->second.size(), MPI_INT, *list_it, MessageType::kCreateGroup1,
                     MPI_COMM_WORLD);
        }
    }
    // 将leader列表发送给所有leader
    std::sort(leader_list.begin(), leader_list.end());
    for (auto it = leader_list.begin(); it != leader_list.end(); ++it) {
        MPI_Send(leader_list.data(), leader_list.size(), MPI_INT, *it, MessageType::kCreateGroup2, MPI_COMM_WORLD);
    }
}

bool Coordinator::CheckWorkerDelay() {
    for (auto it = worker_delay_.begin(); it != worker_delay_.end(); ++it) {
        if (it->second >= max_delay_) {
            LOG(INFO) << "Worker " << it->first << "达到最大延迟，将等待完成更新";
            return false;
        }
    }
    return true;
}

void Coordinator::Terminate() {
    //向所有Worker发送停止命令
    MPI_Status status;
    MPI_Request temp_request;
    for (int i = 0; i < worker_number_; ++i) {
        MPI_Isend(NULL, 0, MPI_INT, i, MessageType::kTerminateCommand, MPI_COMM_WORLD, &temp_request);
        MPI_Request_free(&temp_request);
    }
    //接收到所有Worker的停止报告后，Coordinator也可以退出了
    int count = 0;
    while (count < worker_number_) {
        MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == MessageType::kTerminateComplete) {
            ++count;
        }
    }
    sleep(3);
    LOG(INFO) << "Coordinator停止运行";
}

bool Coordinator::IsSatisfied() {
    if (interval_ > 0 && k_ % interval_ == 0) {
        return ready_worker_list_.size() == worker_number_;
    } else {
        return ready_worker_list_.size() >= min_barrier_ && CheckWorkerDelay();
    }
}

}