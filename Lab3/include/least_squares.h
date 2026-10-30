#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include "utils.h"

#include <string>

void Normal_System(const Vector &x, const Vector &y, int n, Matrix &A, Vector &b);

bool Solve_Gauss(Matrix A, Vector b, Vector &res);

double Poly_Value(const Vector &a, double x);

double Sum_Of_Squares(const Vector &a, const Vector &x, const Vector &y);

void Solve_Task_3_3(const std::string &filename);

#endif // LEAST_SQUARES_H
