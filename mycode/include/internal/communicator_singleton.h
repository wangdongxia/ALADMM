#ifndef SPARSEALLREDUCE_SINGLETON_H
#define SPARSEALLREDUCE_SINGLETON_H

#include <mutex>

namespace spar {

//Communicator单例辅助类
class CommunicatorSingleton {
public:
    static void Set(void *instance);

private:
    static void *instance_;
    static std::mutex instance_mutex_;
};

}

#endif //SPARSEALLREDUCE_SINGLETON_H
