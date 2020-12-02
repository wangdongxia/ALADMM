#ifndef SPARSEALLREDUCE_COORDINATOR_H
#define SPARSEALLREDUCE_COORDINATOR_H

#include <map>
#include <vector>

namespace spar {

class Coordinator {
public:
    static Coordinator *GetInstance();

    void Run();

    void SetMinBarrier(int min_barrier) {
        min_barrier_ = min_barrier;
    }

    void SetMaxDelay(int max_delay) {
        max_delay_ = max_delay;
    }

    void SetInterval(int interval) {
        interval_ = interval;
    }

private:
    int k_;
    int id_;
    int worker_number_;
    int interval_;
    int min_barrier_;
    int max_delay_;
    std::map<int, int> worker_delay_;
    std::vector<int> ready_worker_list_;

    Coordinator();

    void CreateGroup();

    bool CheckWorkerDelay();

    void Terminate();

    bool IsSatisfied();
};

}

#endif //SPARSEALLREDUCE_COORDINATOR_H
