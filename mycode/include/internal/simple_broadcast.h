#ifndef SPARSEALLREDUCE_SIMPLE_BROADCAST_H
#define SPARSEALLREDUCE_SIMPLE_BROADCAST_H

#include <vector>

#include "internal/p2p_communication.h"

namespace spar {

template<class ValueType>
void SimpleBroadcast(ValueType *buffer, int count, int root, int id, std::vector<int> &worker_list, MPI_Comm comm) {
    int worker_number = worker_list.size();
    if (worker_number == 1) return;

    int root_index = -1;
    for (int i = 0; i < worker_number; ++i) {
        if (worker_list[i] == root) {
            root_index = i;
            break;
        }
    }
    CHECK(root_index != -1);

    int sender = root_index;
    int receiver = (root_index + 1) % worker_number;
    for (int i = 0; i < worker_number - 1; ++i) {
        if (id == worker_list[sender]) {
            Send(buffer, count, worker_list[receiver], MessageType::kSimpleBroadcast, comm);
        }
        if (id == worker_list[receiver]) {
            Recv(buffer, count, worker_list[sender], MessageType::kSimpleBroadcast, comm, MPI_STATUS_IGNORE);
        }
        sender = (sender + 1) % worker_number;
        receiver = (receiver + 1) % worker_number;
    }
}

}

#endif //SPARSEALLREDUCE_SIMPLE_BROADCAST_H
