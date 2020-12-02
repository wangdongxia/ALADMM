#include "optimizer/line_search.h"
#include "math/simple_algebra.h"
#include "logging/simple_logging.h"

int BacktrackingLineSearch(Function *function, const double *x, const double *g, const double *d, double *new_x,
                           int dimension, double &step_size, double min_step, double max_step, double c, double r,
                           int max_iterations) {
    double initial_dg = Dot(d, g, dimension);
    if (initial_dg > 0) {
        LOG(WARNING) << "initial d dot g > 0";
        return 1;
    }
    CHECK(step_size > 0);
    if (step_size > max_step) {
        step_size = max_step;
    }
    double function_value, next_function_value;
    function_value = function->Evaluate(x);
    CHECK(!std::isnan(function_value));
    int k = 0;
    while (k < max_iterations) {
        if (step_size < min_step) {
            step_size = min_step;
        }
        for (int i = 0; i < dimension; ++i) {
            new_x[i] = x[i] + step_size * d[i];
        }
        next_function_value = function->Evaluate(new_x);
        CHECK(!std::isnan(next_function_value));
        if (step_size == min_step || next_function_value <= function_value + c * step_size * initial_dg) {
            return 0;
        }
        step_size *= r;
        ++k;
    }
    return 2;
}

