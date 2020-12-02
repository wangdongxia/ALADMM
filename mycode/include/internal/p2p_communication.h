#ifndef SPARSEALLREDUCE_P2P_COMMUNICATION_H
#define SPARSEALLREDUCE_P2P_COMMUNICATION_H

#include <mpi.h>

//以下函数是对MPI函数的封装，采用了模板和模板特化，目的是在编译期间完成相关判定，最大化程序的性能
namespace spar {

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

}

#endif //SPARSEALLREDUCE_P2P_COMMUNICATION_H
