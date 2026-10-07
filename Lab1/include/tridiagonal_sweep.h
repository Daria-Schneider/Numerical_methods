#ifndef TRIDIAGONAL_SWEEP_H
#define TRIDIAGONAL_SWEEP_H

#include <string>
#include "matrix.h"

// a - под главной диагональю (a[0] = 0), b - главная диагональ,
// c - над главной диагональю (c[n-1] = 0), d - правая часть
bool Check_Stability(const Vector &a, const Vector &b, const Vector &c);

bool Tridiagonal_Sweep(const Vector &a, const Vector &b, const Vector &c, const Vector &d,
                       Vector &P, Vector &Q, Vector &x);

void Solve_Task_1_2(const std::string &filename);

#endif // TRIDIAGONAL_SWEEP_H
