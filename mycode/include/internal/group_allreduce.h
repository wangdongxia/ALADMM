#ifndef SPARSEALLREDUCE_GROUP_ALLREDUCE_H
#define SPARSEALLREDUCE_GROUP_ALLREDUCE_H

#include <unistd.h>
#include <vector>

#include "internal/common.h"
#include "internal/reduce_operator.h"
#include "internal/p2p_communication.h"
#include "internal/ring_allreduce.h"
#include "internal/ps_allreduce.h"

namespace spar {

template<class O, class T>
void IntragroupReduce(T *buffer, int count, int id, std::vector<int> &member_list, MPI_Comm comm) {
    //member_list保存了一个组的所有组员id，第一个成员作为这个组的Leader
    int leader_id = member_list.front();
    if (id != leader_id) {
        //普通Worker进程将数据Reduce到Leader进程
        Send(buffer, count, leader_id, MessageType::kIntragroupReduce, comm);
    } else {
        //Leader进程接收所有组员的数据，Reduce到自己的数据上
        MPI_Status status;
        T *recv_buffer = new T[count];
        for (int i = 0; i < member_list.size() - 1; ++i) {
            MPI_Probe(MPI_ANY_SOURCE, MessageType::kIntragroupReduce, comm, &status);
            Recv(recv_buffer, count, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);
            for (int j = 0; j < count; ++j) {
                Reduce<O>(buffer[j], recv_buffer[j]);
            }
        }
        delete[]recv_buffer;
    }
}

template<class T>
void IntragroupBroadcast(T *buffer, int count, int id, std::vector<int> &member_list, MPI_Comm comm) {
    int leader_id = member_list.front();
    if (id != leader_id) {
        //普通Worker在等待Leader广播时，为了不对计算线程产生影响，每睡眠100ms进行一次检查
        int complete = 0;
        MPI_Request request;
        Irecv(buffer, count, leader_id, MessageType::kIntragroupBroadcast, comm, &request);
        MPI_Test(&request, &complete, MPI_STATUS_IGNORE);
        while (!complete) {
            usleep(100000);
            MPI_Test(&request, &complete, MPI_STATUS_IGNORE);
        }
    } else {
        for (auto it = member_list.begin(); it != member_list.end(); ++it) {
            if (*it != id) {
                Send(buffer, count, *it, MessageType::kIntragroupBroadcast, comm);
            }
        }
    }
}

//基于分组的Allreduce算法步骤都相同，分为组内Reduce，组间Allreduce，组内Broadcast，区别在于组间Allreduce的算法不同
template<class O, class T>
void GroupRingAllreduce(T *buffer, int count, int id, std::vector<int> &member_list,
                        std::vector<int> &leader_list, MPI_Comm comm) {
    int leader_id = member_list.front();
    IntragroupReduce<O>(buffer, count, id, member_list, comm);
    if (id == leader_id) {
        //只有Leader进程参与组间Allreduce
        RingAllreduce<O>(buffer, count, id, leader_list, comm);
    }
    IntragroupBroadcast(buffer, count, id, member_list, comm);
}

template<class O, class T>
void SparseGroupRingAllreduce(T *buffer, int count, int id, std::vector<int> &member_list,
                              std::vector<int> &leader_list, MPI_Comm comm) {
    int leader_id = member_list.front();
    IntragroupReduce<O>(buffer, count, id, member_list, comm);
    if (id == leader_id) {
        SparseRingAllreduce<O>(buffer, count, id, leader_list, comm);
    }
    IntragroupBroadcast(buffer, count, id, member_list, comm);
}

template<class O, class T>
void GroupPsAllreduce(T *buffer, int count, int id, std::vector<int> &member_list,
                      std::vector<int> &leader_list, MPI_Comm comm) {
    int leader_id = member_list.front();
    IntragroupReduce<O>(buffer, count, id, member_list, comm);
    if (id == leader_id) {
        PsAllreduce<O>(buffer, count, id, leader_list, comm);
    }
    IntragroupBroadcast(buffer, count, id, member_list, comm);
}

template<class O, class T>
void SparseGroupPsAllreduce(T *buffer, int count, int id, std::vector<int> &member_list,
                            std::vector<int> &leader_list, MPI_Comm comm) {
    int leader_id = member_list.front();
    IntragroupReduce<O>(buffer, count, id, member_list, comm);
    if (id == leader_id) {
        SparsePsAllreduce<O>(buffer, count, id, leader_list, comm);
    }
    IntragroupBroadcast(buffer, count, id, member_list, comm);
}

}

#endif //SPARSEALLREDUCE_GROUP_ALLREDUCE_H
