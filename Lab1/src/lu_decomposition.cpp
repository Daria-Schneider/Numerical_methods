#include "lu_decomposition.h"
#include "file_io.h"

#include <fstream>
#include <utility>

// P * A = L * U, выбор главного элемента по столбцу
bool LU_Decomposition(const Matrix &A, Matrix &P, Matrix &L, Matrix &U, size_t &count_perm)
{
    size_t n = A.size();
    U = A;
    L = Identity_Matrix(n);
    P = Identity_Matrix(n);
    count_perm = 0;

    for (size_t k = 0; k < n; k++)
    {
        // ищем максимальный по модулю элемент в k-м столбце
        size_t n_piv = k;
        for (size_t i = k + 1; i < n; i++)
            if (std::abs(U[i][k]) > std::abs(U[n_piv][k]))
                n_piv = i;

        if (std::abs(U[n_piv][k]) < 1e-12)
            return false; // матрица вырожденная

        if (n_piv != k)
        {
            std::swap(U[k], U[n_piv]);
            std::swap(P[k], P[n_piv]);
            for (size_t j = 0; j < k; j++)
                std::swap(L[k][j], L[n_piv][j]);
            count_perm++;
        }

        // зануляем элементы под диагональю
        for (size_t i = k + 1; i < n; i++)
        {
            double mu = U[i][k] / U[k][k];
            L[i][k] = mu;
            U[i][k] = 0.0;
            for (size_t j = k + 1; j < n; j++)
                U[i][j] -= mu * U[k][j];
        }
    }
    return true;
}

// L * z = b
Vector Forward_Substitution(const Matrix &L, const Vector &b)
{
    size_t n = L.size();
    Vector z(n, 0.0);
    for (size_t i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (size_t j = 0; j < i; j++)
            sum += L[i][j] * z[j];
        z[i] = b[i] - sum;
    }
    return z;
}

// U * x = z
Vector Backward_Substitution(const Matrix &U, const Vector &z)
{
    int n = static_cast<int>(U.size());
    Vector x(n, 0.0);
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
            sum += U[i][j] * x[j];
        x[i] = (z[i] - sum) / U[i][i];
    }
    return x;
}

Vector LU_Solve(const Matrix &P, const Matrix &L, const Matrix &U, const Vector &b)
{
    Vector z = Forward_Substitution(L, Multiply_Vector(P, b));
    return Backward_Substitution(U, z);
}

double LU_Determinant(const Matrix &U, size_t count_perm)
{
    double det = 1.0;
    for (size_t i = 0; i < U.size(); i++)
        det *= U[i][i];

    // каждая перестановка строк меняет знак
    if (count_perm % 2 != 0)
        det = -det;
    return det;
}

// i-й столбец обратной матрицы - решение A * x = e_i
Matrix LU_Inverse(const Matrix &P, const Matrix &L, const Matrix &U)
{
    size_t n = U.size();
    Matrix A_inv(n, Vector(n, 0.0));
    for (size_t i = 0; i < n; i++)
    {
        Vector e(n, 0.0);
        e[i] = 1.0;
        Vector x = LU_Solve(P, L, U, e);
        for (size_t j = 0; j < n; j++)
            A_inv[j][i] = x[j];
    }
    return A_inv;
}

void Solve_Task_1_1(const std::string &filename)
{
    std::cout << "Задание 1.1. LU-разложение\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    size_t n;
    Matrix A;
    Vector b;
    if (!Read_Size(file, n) || !Read_Matrix(file, n, A) || !Read_Vector(file, n, b))
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    std::cout << "Матрица A:\n";
    Print_Matrix(A);
    std::cout << "Вектор b:\n";
    Print_Vector(b);

    Matrix P, L, U;
    size_t count_perm;
    if (!LU_Decomposition(A, P, L, U, count_perm))
    {
        std::cout << "Матрица вырожденная, LU-разложение не существует\n";
        return;
    }

    std::cout << "\nМатрица перестановок P:\n";
    Print_Matrix(P);
    std::cout << "Матрица L:\n";
    Print_Matrix(L);
    std::cout << "Матрица U:\n";
    Print_Matrix(U);

    Matrix LU = Multiply_Matrix(L, U);
    Matrix PA = Multiply_Matrix(P, A);
    std::cout << "L * U:\n";
    Print_Matrix(LU);
    std::cout << "P * A:\n";
    Print_Matrix(PA);
    std::cout << "L * U = P * A: " << (Equal_Matrix(LU, PA) ? "верно" : "неверно") << "\n\n";

    Vector Pb = Multiply_Vector(P, b);
    Vector z = Forward_Substitution(L, Pb);
    Vector x = Backward_Substitution(U, z);
    std::cout << "Прямой ход, z:\n";
    Print_Vector(z);
    std::cout << "Обратный ход, решение x:\n";
    Print_Vector(x);

    std::cout << "\ndet(A) = " << Format_Number(LU_Determinant(U, count_perm)) << "\n\n";

    Matrix A_inv = LU_Inverse(P, L, U);
    std::cout << "Обратная матрица A^(-1):\n";
    Print_Matrix(A_inv);

    Matrix AA_inv = Multiply_Matrix(A, A_inv);
    std::cout << "A * A^(-1):\n";
    Print_Matrix(AA_inv);
    std::cout << "A * A^(-1) = E: " << (Equal_Matrix(AA_inv, Identity_Matrix(n)) ? "верно" : "неверно") << "\n";
}
