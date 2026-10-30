#ifndef NONLINEAR_SYSTEM_H
#define NONLINEAR_SYSTEM_H

#include <string>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// область G: x1_min <= x1 <= x1_max, x2_min <= x2 <= x2_max
struct Region
{
    double x1_min, x1_max;
    double x2_min, x2_max;
};

struct System_Result
{
    Vector x;
    int iterations = 0;
    bool converged = false;
    bool left_region = false;    // приближение вышло из области G
    std::vector<Vector> history; // x(0), x(1), ..., x(k)
    std::vector<double> errors;  // оценка погрешности eps(k) на каждой итерации
};

Vector F_System(const Vector &x, double a);

Matrix Jacobian(const Vector &x, double a);

double Det2(const Matrix &M);

Matrix Inverse2(const Matrix &M);

double Diff_Norm(const Vector &x, const Vector &y);

double Matrix_Norm(const Matrix &M);

bool In_Region(const Vector &x, const Region &G);

Vector Phi_System(const Vector &x, const Matrix &L, double a);

Matrix dPhi_System(const Vector &x, const Matrix &L, double a);

double Find_q_System(const Matrix &L, double a, const Region &G);

System_Result Simple_Iteration_System(const Vector &x0, const Matrix &L, double a, const Region &G, double eps,
                                      double coef);

System_Result Newton_System(const Vector &x0, double a, double eps);

void Solve_Task_2_2(const std::string &filename);

#endif // NONLINEAR_SYSTEM_H
