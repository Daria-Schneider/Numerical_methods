// Задание 3.2. Кубический сплайн с нулевой кривизной на концах: S''(x_0) = S''(x_n) = 0
#include "spline.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Метод прогонки (как в ЛР1, задание 1.2) для трехдиагональной системы
// low_i * c_(i-1) + mid_i * c_i + up_i * c_(i+1) = rhs_i
bool Sweep(const Vector &low, const Vector &mid, const Vector &up, const Vector &rhs, Vector &res)
{
    size_t n = mid.size();
    Vector P(n), Q(n);

    // прямой ход: выражаем c_i = P_i * c_(i+1) + Q_i
    for (size_t i = 0; i < n; i++)
    {
        double denom = mid[i] + (i > 0 ? low[i] * P[i - 1] : 0.0);
        if (std::abs(denom) < 1e-14)
            return false;
        P[i] = -up[i] / denom; // у последней строки up = 0, поэтому P = 0
        Q[i] = (rhs[i] - (i > 0 ? low[i] * Q[i - 1] : 0.0)) / denom;
    }

    // обратный ход
    res.assign(n, 0.0);
    res[n - 1] = Q[n - 1];
    for (int i = static_cast<int>(n) - 2; i >= 0; i--)
        res[i] = P[i] * res[i + 1] + Q[i];
    return true;
}

// Коэффициенты сплайна. Сначала c_2, ..., c_n из системы (3.13),
// потом a_i, b_i, d_i по формулам (3.14).
// system и rhs - матрица и правая часть системы (3.13), их потом выводим на экран
bool Build_Spline(const Vector &x, const Vector &f, Spline &s, Matrix &system, Vector &rhs)
{
    size_t n = x.size() - 1; // число отрезков
    Vector h(n + 1);
    for (size_t i = 1; i <= n; i++)
        h[i] = x[i] - x[i - 1];

    // Система (3.13), уравнение для c_i (i = 2..n) записано в строку i - 2:
    // h_(i-1) c_(i-1) + 2(h_(i-1) + h_i) c_i + h_i c_(i+1) = 3[(f_i - f_(i-1)) / h_i - (f_(i-1) - f_(i-2)) / h_(i-1)]
    // c_1 = 0 (нулевая кривизна слева), поэтому в первом уравнении нет c_1, а в последнем нет c_(n+1).
    // Диагональ 2(h_(i-1) + h_i) больше суммы h_(i-1) + h_i, так что прогонка устойчива
    size_t m = n - 1;
    Vector low(m, 0.0), mid(m, 0.0), up(m, 0.0);
    rhs.assign(m, 0.0);
    for (size_t i = 2; i <= n; i++)
    {
        size_t r = i - 2;
        if (i > 2)
            low[r] = h[i - 1];
        mid[r] = 2 * (h[i - 1] + h[i]);
        if (i < n)
            up[r] = h[i];
        rhs[r] = 3 * ((f[i] - f[i - 1]) / h[i] - (f[i - 1] - f[i - 2]) / h[i - 1]);
    }

    system.assign(m, Vector(m, 0.0));
    for (size_t r = 0; r < m; r++)
    {
        system[r][r] = mid[r];
        if (r > 0)
            system[r][r - 1] = low[r];
        if (r + 1 < m)
            system[r][r + 1] = up[r];
    }

    Vector c_inner;
    if (!Sweep(low, mid, up, rhs, c_inner))
        return false;

    s.x = x;
    s.a.assign(n + 1, 0.0);
    s.b.assign(n + 1, 0.0);
    s.c.assign(n + 2, 0.0); // c_1 = 0 и еще c_(n+1) = 0, про него ниже
    s.d.assign(n + 1, 0.0);
    for (size_t i = 2; i <= n; i++)
        s.c[i] = c_inner[i - 2];

    // Формулы (3.14). Для последнего отрезка в методичке они записаны отдельно:
    // b_n = (f_n - f_(n-1)) / h_n - 2/3 h_n c_n,  d_n = -c_n / (3 h_n).
    // Это те же общие формулы, если взять c_(n+1) = S''(x_n) / 2 = 0, поэтому считаем одним циклом
    for (size_t i = 1; i <= n; i++)
    {
        s.a[i] = f[i - 1];
        s.b[i] = (f[i] - f[i - 1]) / h[i] - h[i] * (s.c[i + 1] + 2 * s.c[i]) / 3;
        s.d[i] = (s.c[i + 1] - s.c[i]) / (3 * h[i]);
    }
    s.c.pop_back(); // c_(n+1) был нужен только для формул
    return true;
}

// номер отрезка i, для которого x_(i-1) <= X <= x_i, или 0, если X вне [x_0, x_n]
size_t Find_Segment(const Spline &s, double X)
{
    for (size_t i = 1; i < s.x.size(); i++)
        if (s.x[i - 1] <= X && X <= s.x[i])
            return i;
    return 0;
}

double Spline_Value(const Spline &s, double X)
{
    size_t i = Find_Segment(s, X);
    double t = X - s.x[i - 1];
    return s.a[i] + s.b[i] * t + s.c[i] * t * t + s.d[i] * t * t * t;
}

// Проверка условий (3.12): S(x_i) = f_i, в каждом внутреннем узле значения S' и S''
// левого и правого куска совпадают, на концах S'' = 0. Возвращает самое большое расхождение
static double Check_Spline(const Spline &s, const Vector &f)
{
    size_t n = s.x.size() - 1;
    double max_diff = std::abs(2 * s.c[1]); // S''(x_0)
    for (size_t i = 1; i <= n; i++)
    {
        // i-й кусок в правом конце отрезка, т.е. в точке x_i
        double h = s.x[i] - s.x[i - 1];
        double S = s.a[i] + s.b[i] * h + s.c[i] * h * h + s.d[i] * h * h * h;
        double dS = s.b[i] + 2 * s.c[i] * h + 3 * s.d[i] * h * h;
        double d2S = 2 * s.c[i] + 6 * s.d[i] * h;

        max_diff = std::max(max_diff, std::abs(s.a[i] - f[i - 1]));
        max_diff = std::max(max_diff, std::abs(S - f[i]));
        if (i < n)
        {
            max_diff = std::max(max_diff, std::abs(dS - s.b[i + 1]));
            max_diff = std::max(max_diff, std::abs(d2S - 2 * s.c[i + 1]));
        }
        else
        {
            max_diff = std::max(max_diff, std::abs(d2S)); // S''(x_n)
        }
    }
    return max_diff;
}

void Solve_Task_3_2(const std::string &filename)
{
    std::cout << "Задание 3.2. Кубический сплайн\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле: точка X*, строка с x_i, строка с f_i
    Vector star, x, f;
    bool ok = Read_Numbers(file, star) && Read_Numbers(file, x) && Read_Numbers(file, f);
    ok = ok && x.size() == f.size() && x.size() >= 3;
    for (size_t i = 1; ok && i < x.size(); i++)
        if (x[i] <= x[i - 1]) // узлы должны идти по возрастанию
            ok = false;
    if (!ok)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }
    double X_star = star[0];
    size_t n = x.size() - 1;

    std::cout << std::setw(5) << "i";
    for (size_t i = 0; i <= n; i++)
        std::cout << std::setw(10) << i;
    std::cout << "\n" << std::setw(5) << "x_i";
    for (size_t i = 0; i <= n; i++)
        std::cout << std::setw(10) << Format_Number(x[i]);
    std::cout << "\n" << std::setw(5) << "f_i";
    for (size_t i = 0; i <= n; i++)
        std::cout << std::setw(10) << Format_Number(f[i]);
    std::cout << "\n\n";

    std::cout << "Шаги h_i = x_i - x_(i-1):";
    for (size_t i = 1; i <= n; i++)
        std::cout << (i > 1 ? ", " : " ") << "h_" << i << " = " << Format_Number(x[i] - x[i - 1]);
    std::cout << "\n";
    std::cout << "Нулевая кривизна на концах: S''(x_0) = 2c_1 = 0, значит c_1 = 0\n\n";

    Spline s;
    Matrix system;
    Vector rhs;
    if (!Build_Spline(x, f, s, system, rhs))
    {
        std::cout << "Прогонка не получилась (деление на ноль)\n";
        return;
    }

    std::vector<std::string> names;
    for (size_t i = 2; i <= n; i++)
        names.push_back("c" + std::to_string(i));

    std::cout << "Система (3.13) для c_2" << (n > 2 ? ", ..., c_" + std::to_string(n) : "") << ":\n";
    for (size_t r = 0; r < system.size(); r++)
        std::cout << "  " << Equation_Text(system[r], names, rhs[r]) << "\n";
    std::cout << "Решение методом прогонки:";
    for (size_t i = 2; i <= n; i++)
        std::cout << (i > 2 ? ", " : " ") << "c_" << i << " = " << Format_Number(s.c[i]);
    std::cout << "\n\n";

    std::cout << "Коэффициенты сплайна (3.14):\n";
    std::cout << std::setw(3) << "i" << std::setw(18) << "[x_(i-1), x_i]" << std::setw(12) << "a_i" << std::setw(12)
              << "b_i" << std::setw(12) << "c_i" << std::setw(12) << "d_i" << "\n";
    for (size_t i = 1; i <= n; i++)
    {
        std::string segment = "[" + Format_Number(x[i - 1]) + ", " + Format_Number(x[i]) + "]";
        std::cout << std::setw(3) << i << std::setw(18) << segment << std::setw(12) << Format_Number(s.a[i])
                  << std::setw(12) << Format_Number(s.b[i]) << std::setw(12) << Format_Number(s.c[i])
                  << std::setw(12) << Format_Number(s.d[i]) << "\n";
    }

    double check = Check_Spline(s, f);
    std::cout << "Проверка: S(x_i) = f_i, S' и S'' непрерывны в узлах, S''(x_0) = S''(x_" << n << ") = 0 - "
              << (check < 1e-9 ? "верно" : "неверно") << "\n\n";

    size_t i = Find_Segment(s, X_star);
    if (i == 0)
    {
        std::cout << "X* = " << X_star << " не лежит на отрезке [" << x[0] << ", " << x[n] << "]\n";
        return;
    }

    std::cout << "X* = " << X_star << " лежит на отрезке [" << x[i - 1] << ", " << x[i] << "], i = " << i << "\n";
    std::string base = Factor_Text(x[i - 1]);
    std::string poly;
    Add_Term(poly, s.a[i], "");
    Add_Term(poly, s.b[i], base);
    Add_Term(poly, s.c[i], base + "^2");
    Add_Term(poly, s.d[i], base + "^3");
    std::cout << "f(x) = " << (poly.empty() ? "0" : poly) << "\n";
    std::cout << "f(" << X_star << ") = " << Format_Number(Spline_Value(s, X_star)) << "\n";
}
