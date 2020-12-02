#ifndef SPARSEALLREDUCE_COLLECTIVE_COMMUNICATION_H
#define SPARSEALLREDUCE_COLLECTIVE_COMMUNICATION_H

#include <mpi.h>
#include <map>

#include "internal/common.h"

namespace spar {

struct ReduceOperator {
};

struct SumOperator : public ReduceOperator {
};

struct MinOperator : public ReduceOperator {
};

struct MaxOperator : public ReduceOperator {
};

struct ProductOperator : public ReduceOperator {
};

template<class T>
inline void Reduce(T &a, const T &b, SumOperator sum) {
    a += b;
}

template<class T>
inline void Reduce(T &a, const T &b, MinOperator min) {
    if (b < a) {
        a = b;
    }
}

template<class T>
inline void Reduce(T &a, const T &b, MaxOperator max) {
    if (b > a) {
        a = b;
    }
}

template<class T>
inline void Reduce(T &a, const T &b, ProductOperator product) {
    a *= b;
}

template<class T>
inline int Send(const T *buf, int count, int dest, int tag, MPI_Comm comm);

template<>
inline int Send(const double *buf, int count, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(buf, count, MPI_DOUBLE, dest, tag, comm);
}

template<>
inline int Send(const float *buf, int count, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(buf, count, MPI_FLOAT, dest, tag, comm);
}

template<>
inline int Send(const int *buf, int count, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(buf, count, MPI_INT, dest, tag, comm);
}

template<>
inline int Send(const char *buf, int count, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(buf, count, MPI_CHAR, dest, tag, comm);
}

template<>
inline int Send(const short *buf, int count, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(buf, count, MPI_SHORT, dest, tag, comm);
}

template<>
inline int Send(const long *buf, int count, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(buf, count, MPI_LONG, dest, tag, comm);
}

template<class T>
inline int Isend(const T *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request);

template<>
inline int Isend(const double *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Isend(buf, count, MPI_DOUBLE, dest, tag, comm, request);
}

template<>
inline int Isend(const float *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Isend(buf, count, MPI_FLOAT, dest, tag, comm, request);
}

template<>
inline int Isend(const int *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Isend(buf, count, MPI_INT, dest, tag, comm, request);
}

template<>
inline int Isend(const char *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Isend(buf, count, MPI_CHAR, dest, tag, comm, request);
}

template<>
inline int Isend(const short *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Isend(buf, count, MPI_SHORT, dest, tag, comm, request);
}

template<>
inline int Isend(const long *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Isend(buf, count, MPI_LONG, dest, tag, comm, request);
}

template<class T>
inline int Recv(T *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status);

template<>
inline int Recv(double *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    return MPI_Recv(buf, count, MPI_DOUBLE, source, tag, comm, status);
}

template<>
inline int Recv(float *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    return MPI_Recv(buf, count, MPI_FLOAT, source, tag, comm, status);
}

template<>
inline int Recv(int *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    return MPI_Recv(buf, count, MPI_INT, source, tag, comm, status);
}

template<>
inline int Recv(char *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    return MPI_Recv(buf, count, MPI_CHAR, source, tag, comm, status);
}

template<>
inline int Recv(short *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    return MPI_Recv(buf, count, MPI_SHORT, source, tag, comm, status);
}

template<>
inline int Recv(long *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    return MPI_Recv(buf, count, MPI_LONG, source, tag, comm, status);
}

template<class T>
inline int Irecv(T *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request);

template<>
inline int Irecv(double *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Irecv(buf, count, MPI_DOUBLE, source, tag, comm, request);
}

template<>
inline int Irecv(float *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Irecv(buf, count, MPI_FLOAT, source, tag, comm, request);
}

template<>
inline int Irecv(int *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Irecv(buf, count, MPI_INT, source, tag, comm, request);
}

template<>
inline int Irecv(char *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Irecv(buf, count, MPI_CHAR, source, tag, comm, request);
}

template<>
inline int Irecv(short *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Irecv(buf, count, MPI_SHORT, source, tag, comm, request);
}

template<>
inline int Irecv(long *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request) {
    return MPI_Irecv(buf, count, MPI_LONG, source, tag, comm, request);
}

template<class OperatorType, class ValueType>
void Allreduce(ValueType *buffer, int count, int id, int process_num, MPI_Comm comm) {
    if (process_num == 1) {
        return;
    }

    OperatorType op;
    ValueType *temp_buf = new ValueType[count];
    if (count < process_num) {
        int sender = 0, receiver = 1;
        for (int i = 0; i < process_num - 1; ++i) {
            if (id == sender) {
                Send(buffer, count, receiver, MessageType::kRingAllreduce, comm);
            }
            if (id == receiver) {
                Recv(temp_buf, count, sender, MessageType::kRingAllreduce, comm, MPI_STATUS_IGNORE);
                for (int j = 0; j < count; ++j) {
                    Reduce(buffer[j], temp_buf[j], op);
                }
            }
            sender = (sender + 1) % process_num;
            receiver = (receiver + 1) % process_num;
        }
        for (int i = 0; i < process_num - 1; ++i) {
            if (id == sender) {
                Send(buffer, count, receiver, MessageType::kRingAllreduce, comm);
            }
            if (id == receiver) {
                Recv(buffer, count, sender, MessageType::kRingAllreduce, comm, MPI_STATUS_IGNORE);
            }
            sender = (sender + 1) % process_num;
            receiver = (receiver + 1) % process_num;
        }
    } else {
        int block_size = count / process_num;
        std::map<int, std::pair<int, int>> blocks;
        for (int i = 0; i < process_num; ++i) {
            blocks[i].first = i * block_size;
            blocks[i].second = block_size;
        }
        blocks[process_num - 1].second = (count - (process_num - 1) * block_size);
        int left = (id - 1 + process_num) % process_num;
        int right = (id + 1) % process_num;
        int send_block_index = id;
        int recv_block_index = (send_block_index - 1 + process_num) % process_num;
        for (int i = 0; i < process_num - 1; ++i) {
            MPI_Request requests[2];
            Isend(buffer + blocks[send_block_index].first, blocks[send_block_index].second, right,
                  MessageType::kRingAllreduceSR, comm, &requests[0]);
            Irecv(temp_buf, blocks[recv_block_index].second, left, MessageType::kRingAllreduceSR, comm, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
            ValueType *base = buffer + blocks[recv_block_index].first;
            block_size = blocks[recv_block_index].second;
            for (int j = 0; j < block_size; ++j) {
                Reduce(base[j], temp_buf[j], op);
            }
            // 下一次循环发送这次循环接收到的segment
            send_block_index = recv_block_index;
            recv_block_index = (recv_block_index - 1 + process_num) % process_num;
        }
        for (int i = 0; i < process_num - 1; ++i) {
            MPI_Request requests[2];
            Isend(buffer + blocks[send_block_index].first, blocks[send_block_index].second, right,
                  MessageType::kRingAllreduceAG, comm, &requests[0]);
            // 由于不需要计算，因此直接覆盖buf上的数据即可
            Irecv(buffer + blocks[recv_block_index].first, blocks[recv_block_index].second, left,
                  MessageType::kRingAllreduceAG, comm, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
            send_block_index = recv_block_index;
            recv_block_index = (recv_block_index - 1 + process_num) % process_num;
        }
    }
    delete[] temp_buf;
}

template<class ValueType>
void Broadcast(ValueType *buffer, int count, int id, int process_num, int root, MPI_Comm comm) {
    if (process_num == 1) {
        return;
    }

    if (count < process_num) {
        int sender = root, receiver = (root + 1) % process_num;
        for (int i = 0; i < process_num - 1; ++i) {
            if (id == sender) {
                Send(buffer, count, receiver, MessageType::kBroadcast, comm);
            }
            if (id == receiver) {
                Recv(buffer, count, sender, MessageType::kBroadcast, comm, MPI_STATUS_IGNORE);
            }
            sender = (sender + 1) % process_num;
            receiver = (receiver + 1) % process_num;
        }
    } else {
        int block_size = count / process_num;
        std::map<int, std::pair<int, int>> blocks;
        for (int i = 0; i < process_num; ++i) {
            blocks[i].first = i * block_size;
            blocks[i].second = block_size;
        }
        blocks[process_num - 1].second = (count - (process_num - 1) * block_size);
        int left = (id - 1 + process_num) % process_num;
        int right = (id + 1) % process_num;
        if (id == root) {
            for (int i = 0; i < process_num - 1; ++i) {
                MPI_Request request;
                Isend(buffer + blocks[i].first, blocks[i].second, right, MessageType::kBroadcast, comm, &request);
                MPI_Request_free(&request);
            }
            Send(buffer + blocks[process_num - 1].first, blocks[process_num - 1].second, right, MessageType::kBroadcast, comm);
        } else {
            for (int i = 0; i < process_num - 1; ++i) {
                Recv(buffer + blocks[i].first, blocks[i].second, left, MessageType::kBroadcast, comm, MPI_STATUS_IGNORE);
                MPI_Request request;
                Isend(buffer + blocks[i].first, blocks[i].second, right, MessageType::kBroadcast, comm, &request);
                MPI_Request_free(&request);
            }
            Recv(buffer + blocks[process_num - 1].first, blocks[process_num - 1].second, left, MessageType::kBroadcast, comm, MPI_STATUS_IGNORE);
            Send(buffer + blocks[process_num - 1].first, blocks[process_num - 1].second, right, MessageType::kBroadcast, comm);
        }
    }
}

}

#endif //SPARSEALLREDUCE_COLLECTIVE_COMMUNICATION_H
