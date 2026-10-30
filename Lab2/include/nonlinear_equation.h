#ifndef NONLINEAR_EQUATION_H
#define NONLINEAR_EQUATION_H

#include <string>
#include <vector>

struct Equation_Result
{
    double x = 0.0;
    int iterations = 0;
    bool converged = false;
    std::vector<double> history; // x(0), x(1), ..., x(k)
    std::vector<double> errors;  // оценка погрешности eps(k) на каждой итерации
};

double F(double x);

double dF(double x);

double d2F(double x);

double Phi(double x);

double dPhi(double x);

double Exact_Root();

double Find_q(double a, double b);

void Phi_Range(double a, double b, double &phi_min, double &phi_max);

bool Keeps_Sign(double (*g)(double), double a, double b);

Equation_Result Simple_Iteration(double x0, double eps, double coef);

Equation_Result Newton_Method(double x0, double eps);

void Solve_Task_2_1(const std::string &filename);

#endif // NONLINEAR_EQUATION_H
