#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "utils.h"

#include <string>

double Y(double x);

double Exact_Integral(double a, double b);

Vector Rectangles(double a, double h, int N);

Vector Trapezoids(double a, double h, int N);

Vector Simpson(double a, double h, int N);

double Runge_Romberg_Error(double F_h, double F_kh, double k, int p);

void Solve_Task_3_5(const std::string &filename);

#endif // INTEGRATION_H
