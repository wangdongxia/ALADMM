#ifndef SPARSEALLREDUCE_SPARSE_ALLREDUCE_H
#define SPARSEALLREDUCE_SPARSE_ALLREDUCE_H

#include <mpi.h>

#include "logging/simple_logging.h"
#include "communicator.h"
#include "coordinator.h"

namespace spar {

inline void Init(int *argc, char ***argv) {
    int provided;
     
    MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
    CHECK_EQ(provided, MPI_THREAD_MULTIPLE) << "当前MPI不支持多线程，请更换MPI版本";
}

inline void Finalize() {
    MPI_Finalize();
}

inline bool IsWorker() {
    int flag;
    MPI_Initialized(&flag);
    CHECK_EQ(flag, 1) << "请先调用spar::Init()";
    int id, process_num;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);
    
    return id != (process_num - 1);
}

}

#endif //SPARSEALLREDUCE_SPARSE_ALLREDUCE_H
