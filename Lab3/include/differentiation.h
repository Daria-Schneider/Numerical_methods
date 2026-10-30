#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H

#include "utils.h"

#include <string>

double First_Derivative_1(const Vector &x, const Vector &y, size_t i);

double First_Derivative_2(const Vector &x, const Vector &y, size_t i, double X);

double Second_Derivative(const Vector &x, const Vector &y, size_t i);

void Solve_Task_3_4(const std::string &filename);

#endif // DIFFERENTIATION_H
