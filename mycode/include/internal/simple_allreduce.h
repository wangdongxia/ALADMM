#ifndef SPARSEALLREDUCE_SIMPLE_ALLREDUCE_H
#define SPARSEALLREDUCE_SIMPLE_ALLREDUCE_H

#include <vector>

#include "internal/reduce_operator.h"
#include "internal/p2p_communication.h"

namespace spar {

//SimpleAllreduce主要在数据量特别少的时候使用
template<class O, class T>
void SimpleAllreduce(T *buffer, int count, int id, std::vector<int> &worker_list, MPI_Comm comm) {
    int worker_number = worker_list.size();
    if (worker_number == 1) return;

    T *recv_buffer = new T[count];
    int sender = 0, receiver = 1;
    for (int i = 0; i < worker_number - 1; ++i) {
        if (id == worker_list[sender]) {
            Send(buffer, count, worker_list[receiver], MessageType::kSimpleAllreduce1, comm);
        }
        if (id == worker_list[receiver]) {
            Recv(recv_buffer, count, worker_list[sender], MessageType::kSimpleAllreduce1, comm, MPI_STATUS_IGNORE);
            for (int j = 0; j < count; ++j) {
                Reduce<O>(buffer[j], recv_buffer[j]);
            }
        }
        sender = (sender + 1) % worker_number;
        receiver = (receiver + 1) % worker_number;
    }

    for (int i = 0; i < worker_number - 1; ++i) {
        if (id == worker_list[sender]) {
            Send(buffer, count, worker_list[receiver], MessageType::kSimpleAllreduce2, comm);
        }
        if (id == worker_list[receiver]) {
            Recv(buffer, count, worker_list[sender], MessageType::kSimpleAllreduce2, comm, MPI_STATUS_IGNORE);
        }
        sender = (sender + 1) % worker_number;
        receiver = (receiver + 1) % worker_number;
    }

    delete []recv_buffer;
}

}

#endif //SPARSEALLREDUCE_SIMPLE_ALLREDUCE_H
