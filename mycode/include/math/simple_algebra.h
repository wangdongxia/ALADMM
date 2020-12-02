#ifndef UTILS_SIMPLE_ALGEBRA_H
#define UTILS_SIMPLE_ALGEBRA_H

#include <cmath>

#include "data/sparse_dataset.h"

void FillZero(double *x, int dimension);

double Dot(const double *x, const double *y, int dimension);

double Dot(const double *x, const Feature *y);

/* L2 Norm */
double Norm(const double *x, int dimension);

double Sigmoid(double x);

/* x = y */
void Assign(double *x, const double *y, int dimension);

/* x = x * scaling_factor */
void Scale(double *x, double scaling_factor, int dimension);

/* x = x + y */
void XPlusY(double *x, const double *y, int dimension);

/* x = ax + y */
void AXPlusY(double *x, double a, const double *y, int dimension);

/* x = x + by */
void XPlusBY(double *x, const double *y, double b, int dimension);

/* x = ax + by */
void AXPlusBY(double *x, double a, const double *y, double b, int dimension);

/* x = x - y */
void XMinusY(double *x, const double *y, int dimension);

/* x = ax - y */
void AXMinusY(double *x, double a, const double *y, int dimension);

/* x = x - by */
void XMinusBY(double *x, const double *y, double b, int dimension);

/* x = ax - by */
void AXMinusBY(double *x, double a, const double *y, double b, int dimension);

#endif
