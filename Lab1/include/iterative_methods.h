#ifndef ITERATIVE_METHODS_H
#define ITERATIVE_METHODS_H

#include <string>
#include "matrix.h"

struct Iteration_Result
{
    Vector x;
    int iterations = 0;
    bool converged = false;
    std::vector<Vector> history; // x(k) на каждой итерации
    std::vector<double> errors;  // оценка погрешности eps(k)
};

bool Has_Diagonal_Dominance(const Matrix &A);

void Build_Equivalent_System(const Matrix &A, const Vector &b, Matrix &alpha, Vector &beta);

int A_Priori_Iterations(double norm_alpha, double norm_beta, double eps);

Iteration_Result Simple_Iterations(const Matrix &alpha, const Vector &beta, double eps, double coef);

Iteration_Result Seidel_Method(const Matrix &alpha, const Vector &beta, double eps, double coef);

void Solve_Task_1_3(const std::string &filename);

#endif // ITERATIVE_METHODS_H
