#ifndef UTILS_LINE_SEARCH_H
#define UTILS_LINE_SEARCH_H

#include "optimizer/differentiable_function.h"

int BacktrackingLineSearch(Function *function, const double *x, const double *g, const double *d, double *new_x,
                           int dimension, double &step_size, double min_step = 1e-15, double max_step = 1e15,
                           double c = 1e-4, double r = 0.5, int max_iterations = 50);

#endif //UTILS_LINE_SEARCH_H
