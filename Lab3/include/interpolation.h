#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "utils.h"

#include <string>

double F(double x);

double dF(double x, int k);

double Omega(const Vector &x, double X);

double dOmega(const Vector &x, size_t i);

double Lagrange(const Vector &x, const Vector &f, double X);

Matrix Divided_Differences(const Vector &x, const Vector &f);

double Newton(const Vector &x, const Matrix &diff, double X);

void Solve_Task_3_1(const std::string &filename);

#endif // INTERPOLATION_H
