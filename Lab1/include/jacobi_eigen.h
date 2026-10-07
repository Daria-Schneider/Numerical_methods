#ifndef JACOBI_EIGEN_H
#define JACOBI_EIGEN_H

#include <string>
#include "matrix.h"

struct Jacobi_Step
{
    int k;
    size_t i, j; // какой элемент зануляли
    double phi;
    double t;    // t(A) после итерации
};

struct Jacobi_Result
{
    Vector eigenvalues;
    Matrix eigenvectors; // по столбцам
    Matrix A_final;
    int iterations = 0;
    bool converged = false;
    double t0 = 0.0;
    std::vector<Jacobi_Step> steps;
};

bool Is_Symmetric(const Matrix &A);

double Off_Diagonal_Norm(const Matrix &A);

Jacobi_Result Jacobi_Rotation(const Matrix &A, double eps);

void Solve_Task_1_4(const std::string &filename);

#endif // JACOBI_EIGEN_H
