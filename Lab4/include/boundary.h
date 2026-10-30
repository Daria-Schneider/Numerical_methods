#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "cauchy.h"
#include "utils.h"

#include <string>
#include <vector>

// одна итерация метода стрельбы: параметр eta, значения y(b), y'(b) и невязка Phi(eta)
struct Shot
{
    double eta, yb, zb, phi;
};

bool Shooting(double h, int N, double eta0, double eta1, double eps, Solution &sol, std::vector<Shot> &shots);

bool Sweep(const Vector &low, const Vector &mid, const Vector &up, const Vector &rhs, Vector &res);

bool Finite_Difference(double h, int N, Vector &x, Vector &y, Vector &low, Vector &mid, Vector &up, Vector &rhs);

void Solve_Task_4_2(const std::string &filename);

#endif // BOUNDARY_H
