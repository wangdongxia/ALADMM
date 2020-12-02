#include "optimizer/line_search.h"
#include "math/simple_algebra.h"
#include "logging/simple_logging.h"
#include "optimizer/lbfgs_optimizer.h"
#include "optimizer/gradient_decent_optimizer.h"

void LbfgsOptimizer::Optimize(double *x) {
    double *d = new double[dimension_];
    double *g = new double[dimension_];
    double *new_g = new double[dimension_];
    double *new_x = new double[dimension_];
    double *alpha = new double[m_];
    double **s = new double *[m_];
    double **y = new double *[m_];
    for (int i = 0; i < m_; ++i) {
        s[i] = new double[dimension_];
        y[i] = new double[dimension_];
    }
    bool failed = false;
    int k = 0, result;
    double step_size, beta, function_value, next_function_value;
    function_->Gradient(x, g);
    function_value = function_->Evaluate(x);
    CHECK(!std::isnan(function_value));
    while (k < max_iterations_) {
        for (int i = 0; i < dimension_; ++i) {
            d[i] = -g[i];
        }
        int index;
        for (int i = k - 1; i >= 0 && i >= k - m_; --i) {
            index = i % m_;
            alpha[index] = Dot(s[index], d, dimension_) / Dot(s[index], y[index], dimension_);
            XMinusBY(d, y[index], alpha[index], dimension_);
        }
        for (int i = (k - m_ > 0 ? k - m_ : 0); i <= k - 1; ++i) {
            index = i % m_;
            beta = Dot(y[index], d, dimension_) / Dot(s[index], y[index], dimension_);
            XPlusBY(d, s[index], alpha[index] - beta, dimension_);
        }
        if (k == 0) {
            step_size = 1;
        } else {
            index = (k - 1) % m_;
            step_size = Dot(s[index], y[index], dimension_) / Dot(y[index], y[index], dimension_);
            step_size = step_size > 0 ? step_size : 1;
        }
        result = BacktrackingLineSearch(function_, x, g, d, new_x, dimension_, step_size);
        if (result == 1) {
            failed = true;
            break;
        } else if (result == 2) {
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
        index = k % m_;
        for (int i = 0; i < dimension_; ++i) {
            y[index][i] = new_g[i] - g[i];
            s[index][i] = new_x[i] - x[i];
        }
        Assign(x, new_x, dimension_);
        Assign(g, new_g, dimension_);
        function_value = next_function_value;
        ++k;
    }
    for (int i = 0; i < m_; ++i) {
        delete[] y[i];
        delete[] s[i];
    }
    delete[] d;
    delete[] g;
    delete[] new_g;
    delete[] new_x;
    delete[] y;
    delete[] s;
    delete[] alpha;

    if (failed) {
        GradientDecentOptimizer optimizer(function_, dimension_, min_gradient_norm_, factor_, max_iterations_);
        optimizer.Optimize(x);
    }
}
