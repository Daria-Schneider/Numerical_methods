#ifndef CAUCHY_H
#define CAUCHY_H

#include "utils.h"

#include <functional>
#include <string>

// Правая часть системы двух уравнений первого порядка (4.15):
// y' = f(x, y, z), z' = g(x, y, z)
using System_Func = std::function<double(double, double, double)>;

// решение задачи Коши в узлах x_k = a + k*h, k = 0..N
struct Solution
{
    Vector x, y, z;
};

Solution Euler(const System_Func &f, const System_Func &g, double a, double y0, double z0, double h, int N);

Solution Runge_Kutta(const System_Func &f, const System_Func &g, double a, double y0, double z0, double h, int N);

Solution Adams(const System_Func &f, const System_Func &g, double a, double y0, double z0, double h, int N);

void Solve_Task_4_1(const std::string &filename);

#endif // CAUCHY_H
