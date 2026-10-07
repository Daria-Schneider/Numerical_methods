#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix Identity_Matrix(size_t n);

Matrix Multiply_Matrix(const Matrix &A, const Matrix &B);

Vector Multiply_Vector(const Matrix &A, const Vector &V);

Matrix Transpose(const Matrix &A);

Vector Subtract_Vector(const Vector &a, const Vector &b);

// max |v_i|
double Vector_Norm(const Vector &v);

// максимальная по строкам сумма |a_ij|
double Matrix_Norm(const Matrix &A);

bool Equal_Matrix(const Matrix &A, const Matrix &B, double eps = 1e-9);

std::string Format_Number(double value, int precision = 6);

std::string Format_Complex(const std::complex<double> &z, int precision = 6);

void Print_Matrix(const Matrix &A);

void Print_Vector(const Vector &V);

#endif // MATRIX_H
