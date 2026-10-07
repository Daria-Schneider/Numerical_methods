#include "tridiagonal_sweep.h"
#include "file_io.h"

#include <fstream>

// достаточные условия устойчивости: a_i, c_i != 0 и |b_i| >= |a_i| + |c_i|,
// причем хотя бы одно неравенство строгое
bool Check_Stability(const Vector &a, const Vector &b, const Vector &c)
{
    size_t n = b.size();
    bool strict = false;
    for (size_t i = 0; i < n; i++)
    {
        if (i > 0 && i + 1 < n && (a[i] == 0 || c[i] == 0))
            return false;

        double side = std::abs(a[i]) + std::abs(c[i]);
        if (std::abs(b[i]) < side)
            return false;
        if (std::abs(b[i]) > side)
            strict = true;
    }
    return strict;
}

bool Tridiagonal_Sweep(const Vector &a, const Vector &b, const Vector &c, const Vector &d,
                       Vector &P, Vector &Q, Vector &x)
{
    size_t n = b.size();
    P.assign(n, 0.0);
    Q.assign(n, 0.0);
    x.assign(n, 0.0);

    // прямой ход
    for (size_t i = 0; i < n; i++)
    {
        double denom = b[i];
        double num = d[i];
        if (i > 0)
        {
            denom += a[i] * P[i - 1];
            num -= a[i] * Q[i - 1];
        }
        if (std::abs(denom) < 1e-14)
            return false;

        P[i] = (i + 1 < n) ? -c[i] / denom : 0.0; // P_n = 0, т.к. c_n = 0
        Q[i] = num / denom;
    }

    // обратный ход
    x[n - 1] = Q[n - 1];
    for (int i = static_cast<int>(n) - 2; i >= 0; i--)
        x[i] = P[i] * x[i + 1] + Q[i];

    return true;
}

void Solve_Task_1_2(const std::string &filename)
{
    std::cout << "Задание 1.2. Метод прогонки\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    size_t n;
    if (!Read_Size(file, n))
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    // храним только три диагонали и правую часть
    Vector a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);
    for (size_t i = 1; i < n; i++)
        file >> a[i];
    for (size_t i = 0; i < n; i++)
        file >> b[i];
    for (size_t i = 0; i + 1 < n; i++)
        file >> c[i];
    for (size_t i = 0; i < n; i++)
        file >> d[i];

    if (!file)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    std::cout << "a (под главной диагональю):\n";
    Print_Vector(a);
    std::cout << "b (главная диагональ):\n";
    Print_Vector(b);
    std::cout << "c (над главной диагональю):\n";
    Print_Vector(c);
    std::cout << "d (правая часть):\n";
    Print_Vector(d);

    if (Check_Stability(a, b, c))
        std::cout << "\nУсловия устойчивости прогонки выполнены\n\n";
    else
        std::cout << "\nУсловия устойчивости прогонки не выполнены, результат может быть неточным\n\n";

    Vector P, Q, x;
    if (!Tridiagonal_Sweep(a, b, c, d, P, Q, x))
    {
        std::cout << "Деление на ноль в прямом ходе, метод неприменим\n";
        return;
    }

    std::cout << "Прогоночные коэффициенты:\n";
    std::cout << std::setw(4) << "i" << std::setw(12) << "P" << std::setw(12) << "Q" << "\n";
    for (size_t i = 0; i < n; i++)
        std::cout << std::setw(4) << i + 1 << std::setw(12) << Format_Number(P[i]) << std::setw(12)
                  << Format_Number(Q[i]) << "\n";

    std::cout << "\nРешение x:\n";
    Print_Vector(x);
}
