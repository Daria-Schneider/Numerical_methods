#ifndef QR_DECOMPOSITION_H
#define QR_DECOMPOSITION_H

#include <complex>
#include <string>
#include "matrix.h"

struct QR_Result
{
    std::vector<std::complex<double>> eigenvalues;
    Matrix A_final;
    int iterations = 0;
    bool converged = false;
};

void QR_Decomposition(const Matrix &A, Matrix &Q, Matrix &R);

QR_Result QR_Algorithm(const Matrix &A, double eps, bool print);

void Solve_Task_1_5(const std::string &filename);

#endif // QR_DECOMPOSITION_H
