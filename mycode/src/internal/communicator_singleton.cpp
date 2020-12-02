#include "internal/communicator_singleton.h"
#include "logging/simple_logging.h"

namespace spar {

std::mutex CommunicatorSingleton::instance_mutex_;

void *CommunicatorSingleton::instance_ = NULL;

void CommunicatorSingleton::Set(void *instance) {
    //只能创建一个Communicator，否则直接结束程序
    if (instance_ == NULL) {
        std::unique_lock<std::mutex> lck(instance_mutex_);
        if (instance_ == NULL) {
            instance_ = instance;
            return;
        }
    }
    CHECK_EQ(instance, instance_) << "只能创建一个Communicator";
}

}
