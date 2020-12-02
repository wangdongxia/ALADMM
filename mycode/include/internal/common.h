#ifndef SPARSEALLREDUCE_COMMON_H
#define SPARSEALLREDUCE_COMMON_H

namespace spar {

//为了避免硬编码，采用enum来定义消息类型
enum MessageType {
    kCreateGroup1 = 0x10000,
    kCreateGroup2,
    kAllreduceApplication,
    kAllreduceCommand,
    kAllreduceComplete,
    kTerminateApplication,
    kTerminateCommand,
    kTerminateComplete,
    kIntragroupReduce,
    kIntragroupBroadcast,
    kSimpleAllreduce1,
    kSimpleAllreduce2,
    kScatterReduce,
    kAllGather,
    kSimpleBroadcast
};

enum class CalculatorState {
    // 已经停止的状态
    kStopped,
    // 希望停止的状态，Communicator会利用这个状态向Coordinator发送停止请求
    kWantToStop,
    // 正在计算
    kCalculating,
    // 计算完成，但还未申请Allreduce
    kCalculationComplete,
    // 申请Allreduce，正在等待结果
    kWaitingAllreduceResult,
    // Allreduce结果已经写入Calculator的缓冲区
    kAllreduceComplete
};

enum class CommunicatorState {
    // 刚初始化的状态
    kWaitToRun,
    // 正在运行的状态
    kRunning,
    // 已经停止的状态
    kStopped
};

}

#endif //SPARSEALLREDUCE_COMMON_H
