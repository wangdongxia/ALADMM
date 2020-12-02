#ifndef UTILS_DIFFERENTIABLE_FUNCTION_H
#define UTILS_DIFFERENTIABLE_FUNCTION_H

class Function {
public:
    virtual double Evaluate(const double *x) = 0;
};

class DifferentiableFunction : public Function {
public:
    virtual void Gradient(const double *x, double *g) = 0;
};

#endif //DIFFERENTIABLE_FUNCTION_H
