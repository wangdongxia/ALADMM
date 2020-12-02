#ifndef SPARSEALLREDUCE_PS_ALLREDUCE_H
#define SPARSEALLREDUCE_PS_ALLREDUCE_H

#include <map>
#include <vector>
#include <cstring>

#include "logging/simple_logging.h"
#include "internal/common.h"
#include "internal/reduce_operator.h"
#include "internal/simple_allreduce.h"
#include "internal/p2p_communication.h"

//PS-Allreduce和稀疏PS-Allreduce的算法实现
namespace spar {

template<class O, class T>
void PsAllreduce(T *buffer, int count, int id, std::vector<int> &worker_list, MPI_Comm comm) {
    //worker_list包含了参与此次Allreduce的所有进程id
    int worker_number = worker_list.size();
    if (worker_number == 1) return;
    //如果数据量特别少的话，直接使用SimpleAllreduce
    if (count < 512 * worker_number) {
        SimpleAllreduce<O>(buffer, count, id, worker_list, comm);
        return;
    }

    //对数组进行分块，在blocks中记录每一块的起始地址和块大小
    int block_size = count / worker_number;
    std::map<int, std::pair<int, int>> blocks;
    for (int i = 0; i < worker_number; ++i) {
        blocks[i].first = i * block_size;
        blocks[i].second = block_size;
    }
    blocks[worker_number - 1].second = (count - (worker_number - 1) * block_size);
    //确定本进程在worker_list中的索引，从而确定了本进程管理的块
    int my_index = -1;
    for (int i = 0; i < worker_number; ++i) {
        if (worker_list[i] == id) {
            my_index = i;
            break;
        }
    }
    CHECK_NE(my_index, -1);

    MPI_Status status;
    T *base;
    T *recv_buffer = new T[count];
    std::vector<MPI_Request> request_list;
    for (int i = 0; i < worker_number; ++i) {
        //除了自己管理的块以外，将其他块发送到对应的进程，采用非阻塞发送
        if (i != my_index) {
            request_list.emplace_back();
            base = buffer + blocks[i].first;
            block_size = blocks[i].second;
            Isend(base, block_size, worker_list[i], MessageType::kScatterReduce, comm, &request_list.back());
        }
    }
    base = buffer + blocks[my_index].first;
    block_size = blocks[my_index].second;
    for (int i = 0; i < worker_number - 1; ++i) {
        //从其他进程接收本进程管理的块，并Reduce到本进程的数组上
        MPI_Probe(MPI_ANY_SOURCE, MessageType::kScatterReduce, comm, &status);
        Recv(recv_buffer, block_size, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);
        for (int j = 0; j < block_size; ++j) {
            Reduce<O>(base[j], recv_buffer[j]);
        }
    }
    MPI_Waitall(request_list.size(), request_list.data(), MPI_STATUSES_IGNORE);

    request_list.clear();
    for (int i = 0; i < worker_number; ++i) {
        //把自己管理的块发送到其他进程
        if (i != my_index) {
            request_list.emplace_back();
            Isend(base, block_size, worker_list[i], MessageType::kAllGather, comm, &request_list.back());
        }
    }
    for (int i = 0; i < worker_number - 1; ++i) {
        //接收其他进程管理的块，并覆盖本进程上的数据
        MPI_Probe(MPI_ANY_SOURCE, MessageType::kAllGather, comm, &status);
        int sender_index = -1;
        for (int j = 0; j < worker_number; ++j) {
            if (worker_list[j] == status.MPI_SOURCE) {
                sender_index = j;
                break;
            }
        }
        CHECK_NE(sender_index, -1);
        base = buffer + blocks[sender_index].first;
        block_size = blocks[sender_index].second;
        Recv(base, block_size, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);
    }
    MPI_Waitall(request_list.size(), request_list.data(), MPI_STATUSES_IGNORE);

    delete[] recv_buffer;
}

//稀疏版本的算法和稠密版本的步骤基本一致，只是在发送前判断是否满足稀疏传输条件
template<class O, class T>
void SparsePsAllreduce(T *buffer, int count, int id, std::vector<int> &worker_list, MPI_Comm comm) {
    int worker_number = worker_list.size();
    if (worker_number == 1) return;
    if (count < 512 * worker_number) {
        SimpleAllreduce<O>(buffer, count, id, worker_list, comm);
        return;
    }

    int block_size = count / worker_number;
    std::map<int, std::pair<int, int>> blocks;
    for (int i = 0; i < worker_number; ++i) {
        blocks[i].first = i * block_size;
        blocks[i].second = block_size;
    }
    blocks[worker_number - 1].second = (count - (worker_number - 1) * block_size);
    int my_index = -1;
    for (int i = 0; i < worker_number; ++i) {
        if (worker_list[i] == id) {
            my_index = i;
            break;
        }
    }
    CHECK_NE(my_index, -1);

    T *base;
    MPI_Status status;
    std::vector<MPI_Request> request_list;
    bool need_check = true;
    int isize = sizeof(int);
    int vsize = sizeof(T);
    int *index_buffer = new int[count];
    T *value_buffer = new T[count];
    T *recv_buffer = new T[count];
    for (int i = 0; i < worker_number; ++i) {
        if (i != my_index) {
            request_list.emplace_back();
            base = buffer + blocks[i].first;
            block_size = blocks[i].second;
            int nnz = 0;
            for (int j = 0; j < block_size; ++j) {
                if (base[j] != 0) {
                    ++nnz;
                }
            }
            //如果满足稀疏传输条件则采用稀疏传输
            if (nnz * (isize + vsize) < block_size * vsize) {
                int k = 0;
                T *base2 = value_buffer + blocks[i].first;
                for (int j = 0; j < block_size; ++j) {
                    if (base[j] != 0) {
                        base2[k] = base[j];
                        index_buffer[k] = j;
                        ++k;
                    }
                }
                //把元素的值和元素的索引放到同一个内存空间上，一同发送，不分两次发送了
                memcpy(base2 + nnz, index_buffer, isize * nnz);
                MPI_Isend(base2, nnz * (isize + vsize), MPI_CHAR, worker_list[i],
                          MessageType::kScatterReduce, comm, &request_list.back());
            } else {
                MPI_Isend(base, block_size * vsize, MPI_CHAR, worker_list[i], MessageType::kScatterReduce, comm,
                          &request_list.back());
            }
        }
    }
    base = buffer + blocks[my_index].first;
    block_size = blocks[my_index].second;
    for (int i = 0; i < worker_number - 1; ++i) {
        MPI_Recv(recv_buffer, block_size * vsize, MPI_CHAR, MPI_ANY_SOURCE, MessageType::kScatterReduce, comm, &status);
        int nnz = 0;
        MPI_Get_count(&status, MPI_CHAR, &nnz);
        //如果接收到的字节数少于这个块应有的大小，说明进程发送的是稀疏数据
        if (nnz < block_size * vsize) {
            CHECK_EQ(nnz % (isize + vsize), 0);
            nnz = nnz / (isize + vsize);
            memcpy(index_buffer, recv_buffer + nnz, nnz * isize);
            for (int j = 0; j < nnz; ++j) {
                Reduce<O>(base[index_buffer[j]], recv_buffer[j]);
            }
        } else {
            CHECK_EQ(nnz, block_size * vsize);
            for (int j = 0; j < block_size; ++j) {
                Reduce<O>(base[j], recv_buffer[j]);
            }
            //如果接收到稠密的块说明之后发送本进程管理的块时不再需要检查是否满足稀疏传输条件，必然稠密
            need_check = false;
        }
    }
    MPI_Waitall(request_list.size(), request_list.data(), MPI_STATUSES_IGNORE);

    request_list.clear();
    if (need_check) {
        int nnz = 0;
        for (int j = 0; j < block_size; ++j) {
            if (base[j] != 0) {
                ++nnz;
            }
        }
        if (nnz * (isize + vsize) < block_size * vsize) {
            int k = 0;
            for (int j = 0; j < block_size; ++j) {
                if (base[j] != 0) {
                    value_buffer[k] = base[j];
                    index_buffer[k] = j;
                    ++k;
                }
            }
            memcpy(value_buffer + nnz, index_buffer, isize * nnz);
            for (int i = 0; i < worker_number; ++i) {
                if (i != my_index) {
                    request_list.emplace_back();
                    MPI_Isend(value_buffer, nnz * (isize + vsize), MPI_CHAR, worker_list[i], MessageType::kAllGather,
                              comm, &request_list.back());
                }
            }
        } else {
            need_check = false;
        }
    }
    if (!need_check) {
        for (int i = 0; i < worker_number; ++i) {
            if (i != my_index) {
                request_list.emplace_back();
                MPI_Isend(base, block_size * vsize, MPI_CHAR, worker_list[i], MessageType::kAllGather, comm,
                          &request_list.back());
            }
        }
    }
    for (int i = 0; i < worker_number - 1; ++i) {
        MPI_Probe(MPI_ANY_SOURCE, MessageType::kAllGather, comm, &status);
        int sender_index = -1;
        for (int j = 0; j < worker_number; ++j) {
            if (worker_list[j] == status.MPI_SOURCE) {
                sender_index = j;
                break;
            }
        }
        CHECK_NE(sender_index, -1);

        base = buffer + blocks[sender_index].first;
        block_size = blocks[sender_index].second;
        MPI_Recv(recv_buffer, block_size * vsize, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm,
                 MPI_STATUS_IGNORE);
        int nnz = 0;
        MPI_Get_count(&status, MPI_CHAR, &nnz);
        if (nnz < block_size * vsize) {
            CHECK_EQ(nnz % (isize + vsize), 0);
            nnz = nnz / (isize + vsize);
            memcpy(index_buffer, recv_buffer + nnz, nnz * isize);
            for (int j = 0; j < block_size; ++j) {
                base[j] = 0;
            }
            for (int j = 0; j < nnz; ++j) {
                base[index_buffer[j]] = recv_buffer[j];
            }
        } else {
            CHECK_EQ(nnz, block_size * vsize);
            for (int j = 0; j < block_size; ++j) {
                base[j] = recv_buffer[j];
            }
        }
    }
    MPI_Waitall(request_list.size(), request_list.data(), MPI_STATUSES_IGNORE);

    delete[] recv_buffer;
    delete[] value_buffer;
    delete[] index_buffer;
}

}

#endif //SPARSEALLREDUCE_PS_ALLREDUCE_H
