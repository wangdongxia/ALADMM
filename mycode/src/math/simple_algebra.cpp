#include "math/simple_algebra.h"

void FillZero(double *x, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] = 0;
    }
}

double Dot(const double *x, const double *y, int dimension) {
    double sum = 0;
    for (int i = 0; i < dimension; ++i) {
        sum += x[i] * y[i];
    }
    return sum;
}

double Dot(const double *x, const Feature *y) {
    double sum = 0;
    while (y->index != -1) {
        sum += x[y->index] * y->value;
        ++y;
    }
    return sum;
}

/* L2 Norm */
double Norm(const double *x, int dimension) {
    return std::sqrt(Dot(x, x, dimension));
}

double Sigmoid(const double x) {
    return 1.0 / (1 + std::exp(-x));
}

/* x = y */
void Assign(double *x, const double *y, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] = y[i];
    }
}

/* x = x * scaling_factor */
void Scale(double *x, double scaling_factor, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] *= scaling_factor;
    }
}

/* x = x + y */
void XPlusY(double *x, const double *y, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] += y[i];
    }
}

/* x = ax + y */
void AXPlusY(double *x, double a, const double *y, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] = a * x[i] + y[i];
    }
}

/* x = x + by */
void XPlusBY(double *x, const double *y, double b, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] += b * y[i];
    }
}

/* x = ax + by */
void AXPlusBY(double *x, double a, const double *y, double b, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] = a * x[i] + b * y[i];
    }
}

/* x = x - y */
void XMinusY(double *x, const double *y, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] -= y[i];
    }
}

/* x = ax - y */
void AXMinusY(double *x, double a, const double *y, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] = a * x[i] - y[i];
    }
}

/* x = x - by */
void XMinusBY(double *x, const double *y, double b, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] -= b * y[i];
    }
}

/* x = ax - by */
void AXMinusBY(double *x, double a, const double *y, double b, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        x[i] = a * x[i] - b * y[i];
    }
}
