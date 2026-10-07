#include "matrix.h"

#include <algorithm>

Matrix Identity_Matrix(size_t n)
{
    Matrix E(n, Vector(n, 0.0));
    for (size_t i = 0; i < n; i++)
        E[i][i] = 1.0;
    return E;
}

Matrix Multiply_Matrix(const Matrix &A, const Matrix &B)
{
    if (A.empty() || B.empty())
        return {};

    size_t rows = A.size();
    size_t cols = B[0].size();
    size_t inner = A[0].size();

    Matrix res(rows, Vector(cols, 0.0));
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            for (size_t k = 0; k < inner; k++)
                res[i][j] += A[i][k] * B[k][j];
    return res;
}

Vector Multiply_Vector(const Matrix &A, const Vector &V)
{
    Vector res(A.size(), 0.0);
    for (size_t i = 0; i < A.size(); i++)
        for (size_t k = 0; k < V.size(); k++)
            res[i] += A[i][k] * V[k];
    return res;
}

Matrix Transpose(const Matrix &A)
{
    if (A.empty())
        return {};

    Matrix T(A[0].size(), Vector(A.size(), 0.0));
    for (size_t i = 0; i < A.size(); i++)
        for (size_t j = 0; j < A[0].size(); j++)
            T[j][i] = A[i][j];
    return T;
}

Vector Subtract_Vector(const Vector &a, const Vector &b)
{
    Vector res(a.size());
    for (size_t i = 0; i < a.size(); i++)
        res[i] = a[i] - b[i];
    return res;
}

double Vector_Norm(const Vector &v)
{
    double norm = 0.0;
    for (double value : v)
        norm = std::max(norm, std::abs(value));
    return norm;
}

double Matrix_Norm(const Matrix &A)
{
    double norm = 0.0;
    for (const auto &row : A)
    {
        double sum = 0.0;
        for (double value : row)
            sum += std::abs(value);
        norm = std::max(norm, sum);
    }
    return norm;
}

bool Equal_Matrix(const Matrix &A, const Matrix &B, double eps)
{
    if (A.size() != B.size())
        return false;

    for (size_t i = 0; i < A.size(); i++)
    {
        if (A[i].size() != B[i].size())
            return false;
        for (size_t j = 0; j < A[i].size(); j++)
            if (std::abs(A[i][j] - B[i][j]) > eps)
                return false;
    }
    return true;
}

std::string Format_Number(double value, int precision)
{
    // погрешности округления вида 1e-17 печатаем как 0
    if (std::abs(value) < 1e-10)
        value = 0.0;

    std::ostringstream oss;
    oss << std::setprecision(precision) << value;
    return oss.str();
}

std::string Format_Complex(const std::complex<double> &z, int precision)
{
    std::string res = Format_Number(z.real(), precision);
    if (z.imag() != 0.0)
    {
        res += (z.imag() > 0 ? " + " : " - ");
        res += Format_Number(std::abs(z.imag()), precision) + "i";
    }
    return res;
}

void Print_Matrix(const Matrix &A)
{
    for (const auto &row : A)
    {
        for (double value : row)
            std::cout << std::setw(13) << Format_Number(value);
        std::cout << '\n';
    }
    std::cout << '\n';
}

void Print_Vector(const Vector &V)
{
    for (double value : V)
        std::cout << std::setw(13) << Format_Number(value);
    std::cout << '\n';
}
