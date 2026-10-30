// Задание 3.3. Метод наименьших квадратов, приближающие многочлены 1-й и 2-й степени
#include "least_squares.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Нормальная система МНК (3.17) для многочлена F_n(x) = a_0 + a_1 x + ... + a_n x^n:
// sum_i a_i * sum_j x_j^(k+i) = sum_j y_j * x_j^k, k = 0..n.
// Она получается из условия минимума суммы квадратов (3.15): все производные dФ/da_k = 0
void Normal_System(const Vector &x, const Vector &y, int n, Matrix &A, Vector &b)
{
    A.assign(n + 1, Vector(n + 1, 0.0));
    b.assign(n + 1, 0.0);
    for (int k = 0; k <= n; k++)
    {
        for (int i = 0; i <= n; i++)
            for (size_t j = 0; j < x.size(); j++)
                A[k][i] += std::pow(x[j], k + i);
        for (size_t j = 0; j < x.size(); j++)
            b[k] += y[j] * std::pow(x[j], k);
    }
}

// Метод Гаусса с выбором главного элемента по столбцу
bool Solve_Gauss(Matrix A, Vector b, Vector &res)
{
    size_t n = b.size();
    for (size_t col = 0; col < n; col++)
    {
        // ищем строку с самым большим по модулю элементом в столбце и ставим ее наверх
        size_t best = col;
        for (size_t row = col + 1; row < n; row++)
            if (std::abs(A[row][col]) > std::abs(A[best][col]))
                best = row;
        if (std::abs(A[best][col]) < 1e-14)
            return false;
        std::swap(A[col], A[best]);
        std::swap(b[col], b[best]);

        // зануляем столбец под главным элементом
        for (size_t row = col + 1; row < n; row++)
        {
            double m = A[row][col] / A[col][col];
            for (size_t k = col; k < n; k++)
                A[row][k] -= m * A[col][k];
            b[row] -= m * b[col];
        }
    }

    // обратный ход
    res.assign(n, 0.0);
    for (int i = static_cast<int>(n) - 1; i >= 0; i--)
    {
        double sum = b[i];
        for (size_t k = i + 1; k < n; k++)
            sum -= A[i][k] * res[k];
        res[i] = sum / A[i][i];
    }
    return true;
}

// значение многочлена a_0 + a_1 x + ... + a_n x^n по схеме Горнера
double Poly_Value(const Vector &a, double x)
{
    double value = 0.0;
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; i--)
        value = value * x + a[i];
    return value;
}

// сумма квадратов ошибок Ф = sum (F(x_j) - y_j)^2, формула (3.15)
double Sum_Of_Squares(const Vector &a, const Vector &x, const Vector &y)
{
    double sum = 0.0;
    for (size_t j = 0; j < x.size(); j++)
    {
        double e = Poly_Value(a, x[j]) - y[j];
        sum += e * e;
    }
    return sum;
}

void Solve_Task_3_3(const std::string &filename)
{
    std::cout << "Задание 3.3. Метод наименьших квадратов\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле две строки: x_j и y_j
    Vector x, y;
    if (!Read_Numbers(file, x) || !Read_Numbers(file, y) || x.size() != y.size() || x.size() < 3)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }
    size_t N = x.size() - 1; // точки нумеруются от 0 до N, как в методичке

    std::cout << std::setw(5) << "i";
    for (size_t j = 0; j <= N; j++)
        std::cout << std::setw(10) << j;
    std::cout << "\n" << std::setw(5) << "x_i";
    for (size_t j = 0; j <= N; j++)
        std::cout << std::setw(10) << Format_Number(x[j]);
    std::cout << "\n" << std::setw(5) << "y_i";
    for (size_t j = 0; j <= N; j++)
        std::cout << std::setw(10) << Format_Number(y[j]);
    std::cout << "\n\n";

    Vector coefs[3]; // coefs[n] - коэффициенты многочлена степени n
    double phi[3] = {0, 0, 0};
    for (int n = 1; n <= 2; n++)
    {
        std::cout << "Многочлен " << n << "-й степени: F" << n << "(x) = a0 + a1*x" << (n == 2 ? " + a2*x^2" : "")
                  << "\n";

        Matrix A;
        Vector b;
        Normal_System(x, y, n, A, b);

        std::vector<std::string> names;
        for (int i = 0; i <= n; i++)
            names.push_back("a" + std::to_string(i));

        std::cout << "Нормальная система (3.17):\n";
        for (int k = 0; k <= n; k++)
            std::cout << "  " << Equation_Text(A[k], names, b[k]) << "\n";

        if (!Solve_Gauss(A, b, coefs[n]))
        {
            std::cout << "Система вырожденная, решить нельзя\n";
            return;
        }

        std::cout << "Решение методом Гаусса:";
        for (int i = 0; i <= n; i++)
            std::cout << (i > 0 ? ", " : " ") << "a" << i << " = " << Format_Number(coefs[n][i]);
        std::cout << "\n";

        std::string poly;
        Add_Term(poly, coefs[n][0], "");
        Add_Term(poly, coefs[n][1], "x");
        if (n == 2)
            Add_Term(poly, coefs[n][2], "x^2");
        std::cout << "F" << n << "(x) = " << (poly.empty() ? "0" : poly) << "\n\n";

        std::cout << std::setw(3) << "i" << std::setw(10) << "x_i" << std::setw(12) << "y_i" << std::setw(12)
                  << ("F" + std::to_string(n) + "(x_i)") << std::setw(14) << "F(x_i) - y_i" << "\n";
        for (size_t j = 0; j <= N; j++)
        {
            double F = Poly_Value(coefs[n], x[j]);
            std::cout << std::setw(3) << j << std::setw(10) << Format_Number(x[j]) << std::setw(12)
                      << Format_Number(y[j]) << std::setw(12) << Format_Number(F) << std::setw(14)
                      << Format_Number(F - y[j]) << "\n";
        }

        phi[n] = Sum_Of_Squares(coefs[n], x, y);
        std::cout << "Сумма квадратов ошибок Ф" << n << " = " << Format_Number(phi[n]) << "\n\n";
    }

    std::cout << "Ф1 = " << Format_Number(phi[1]) << ", Ф2 = " << Format_Number(phi[2]) << "\n";
    if (phi[2] < phi[1])
        std::cout << "Многочлен 2-й степени приближает таблицу лучше\n";
    else
        std::cout << "Многочлен 2-й степени приближает таблицу не лучше, чем 1-й\n";
    std::cout << "Графики строит скрипт Lab3/plots.py\n";
}
