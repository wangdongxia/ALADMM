#include "optimizer/line_search.h"
#include "math/simple_algebra.h"
#include "logging/simple_logging.h"
#include "optimizer/gradient_decent_optimizer.h"

void GradientDecentOptimizer::Optimize(double *x) {
    double *d = new double[dimension_];
    double *g = new double[dimension_];
    double *new_g = new double[dimension_];
    double *new_x = new double[dimension_];
    double *s = new double[dimension_];
    double *y = new double[dimension_];

    int k = 0, result;
    double step_size, function_value, next_function_value;
    function_->Gradient(x, g);
    function_value = function_->Evaluate(x);
    CHECK(!std::isnan(function_value));
    while (k < max_iterations_) {
        for (int i = 0; i < dimension_; ++i) {
            d[i] = -g[i];
        }
        if (k == 0) {
            step_size = 1;
        } else {
            step_size = Dot(s, y, dimension_) / Dot(y, y, dimension_);
            step_size = step_size > 0 ? step_size : 1;
        }
        result = BacktrackingLineSearch(function_, x, g, d, new_x, dimension_, step_size);
        if (result != 0) {
            for (int i = 0; i < dimension_; ++i) {
                new_x[i] = x[i] + 1e-6 * d[i];
            }
        }
        function_->Gradient(new_x, new_g);
        next_function_value = function_->Evaluate(new_x);
        CHECK(!std::isnan(next_function_value));
        if (Norm(new_g, dimension_) <= min_gradient_norm_) {
            break;
        }
        double denom = std::max(std::max(std::abs(function_value), std::abs(next_function_value)), 1.0);
        if ((function_value - next_function_value) / denom <= factor_) {
            break;
        }
        for (int i = 0; i < dimension_; ++i) {
            s[i] = new_x[i] - x[i];
            y[i] = new_g[i] - g[i];
        }
        Assign(x, new_x, dimension_);
        Assign(g, new_g, dimension_);
        function_value = next_function_value;
        ++k;
    }

    delete[] d;
    delete[] g;
    delete[] new_g;
    delete[] new_x;
    delete[] y;
    delete[] s;
}

