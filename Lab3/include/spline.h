#ifndef SPLINE_H
#define SPLINE_H

#include "utils.h"

#include <string>

// Кубический сплайн (3.11): на отрезке [x_(i-1), x_i]
// S(x) = a_i + b_i (x - x_(i-1)) + c_i (x - x_(i-1))^2 + d_i (x - x_(i-1))^3, i = 1..n.
// Нумерация как в методичке, с единицы, поэтому элемент [0] у a, b, c, d не используется
struct Spline
{
    Vector x;
    Vector a, b, c, d;
};

bool Sweep(const Vector &low, const Vector &mid, const Vector &up, const Vector &rhs, Vector &res);

bool Build_Spline(const Vector &x, const Vector &f, Spline &s, Matrix &system, Vector &rhs);

size_t Find_Segment(const Spline &s, double X);

double Spline_Value(const Spline &s, double X);

void Solve_Task_3_2(const std::string &filename);

#endif // SPLINE_H
