// Задание 4.1. Задача Коши для ОДУ 2-го порядка: методы Эйлера, Рунге-Кутты и Адамса 4-го порядка,
// оценка погрешности методом Рунге-Ромберга и сравнением с точным решением
#include "cauchy.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

const double PI = std::acos(-1.0);

// Вариант 27: x^2 y'' - 2x y' + (x^2 + 2) y = 0, y(pi/2) = pi/2, y'(pi/2) = 1 - pi/2, x in [pi/2, pi/2 + 1].
// Все методы написаны для уравнений первого порядка, поэтому вводим z = y' и получаем систему (4.21):
//   y' = z
//   z' = (2x z - (x^2 + 2) y) / x^2   (это y'', выраженная из уравнения; на нашем отрезке x != 0)
// с начальными условиями y(pi/2) = pi/2, z(pi/2) = 1 - pi/2
static double F(double, double, double z)
{
    return z;
}

static double G(double x, double y, double z)
{
    return (2 * x * z - (x * x + 2) * y) / (x * x);
}

// точное решение, только для проверки
static double Exact(double x)
{
    return x * std::cos(x) + x * std::sin(x);
}

// сетка x_k = a + k*h, k = 0..N, и начальные условия в узле x_0
static Solution Make_Grid(double a, double y0, double z0, double h, int N)
{
    Solution s;
    s.x.resize(N + 1);
    for (int k = 0; k <= N; k++)
        s.x[k] = a + k * h;
    s.y.assign(N + 1, 0.0);
    s.z.assign(N + 1, 0.0);
    s.y[0] = y0;
    s.z[0] = z0;
    return s;
}

// Явный метод Эйлера (4.2): y_(k+1) = y_k + h f(x_k, y_k).
// Идем по касательной: производная в начале шага известна из уравнения. Для системы делаем
// это для y и для z одновременно. Первый порядок точности: при уменьшении h в 2 раза
// погрешность уменьшается тоже примерно в 2 раза
Solution Euler(const System_Func &f, const System_Func &g, double a, double y0, double z0, double h, int N)
{
    Solution s = Make_Grid(a, y0, z0, h, N);
    for (int k = 0; k < N; k++)
    {
        s.y[k + 1] = s.y[k] + h * f(s.x[k], s.y[k], s.z[k]);
        s.z[k + 1] = s.z[k] + h * g(s.x[k], s.y[k], s.z[k]);
    }
    return s;
}

// Метод Рунге-Кутты 4-го порядка для системы (4.16).
// На каждом шаге считаем наклон в четырех точках: в начале шага (K1, L1), дважды в середине (K2, L2 и K3, L3)
// и в конце (K4, L4). K - приращения для y, L - для z. Берем их среднее с весами 1, 2, 2, 1.
// Это совпадает с рядом Тейлора до h^4, поэтому порядок точности 4
Solution Runge_Kutta(const System_Func &f, const System_Func &g, double a, double y0, double z0, double h, int N)
{
    Solution s = Make_Grid(a, y0, z0, h, N);
    for (int k = 0; k < N; k++)
    {
        double x = s.x[k], y = s.y[k], z = s.z[k];

        double K1 = h * f(x, y, z);
        double L1 = h * g(x, y, z);
        double K2 = h * f(x + h / 2, y + K1 / 2, z + L1 / 2);
        double L2 = h * g(x + h / 2, y + K1 / 2, z + L1 / 2);
        double K3 = h * f(x + h / 2, y + K2 / 2, z + L2 / 2);
        double L3 = h * g(x + h / 2, y + K2 / 2, z + L2 / 2);
        double K4 = h * f(x + h, y + K3, z + L3);
        double L4 = h * g(x + h, y + K3, z + L3);

        s.y[k + 1] = y + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
        s.z[k + 1] = z + (L1 + 2 * L2 + 2 * L3 + L4) / 6;
    }
    return s;
}

// Метод Адамса 4-го порядка (4.25):
//   y_(k+1) = y_k + h/24 * (55 f_k - 59 f_(k-1) + 37 f_(k-2) - 9 f_(k-3)),  f_i = f(x_i, y_i, z_i)
// Откуда формула: y_(k+1) = y_k + интеграл от x_k до x_(k+1) от f(x, y(x)), и f под интегралом
// заменяем многочленом 3-й степени, построенным по четырем последним узлам.
// Метод многошаговый и сам стартовать не может: для первого шага уже нужны y_0, y_1, y_2, y_3.
// Поэтому y_1, y_2, y_3 считаем методом Рунге-Кутты того же 4-го порядка, чтобы не потерять точность.
// Плюс метода: на шаге правая часть вычисляется только один раз (в Рунге-Кутте четыре раза)
Solution Adams(const System_Func &f, const System_Func &g, double a, double y0, double z0, double h, int N)
{
    Solution s = Make_Grid(a, y0, z0, h, N);

    int start = std::min(N, 3);
    Solution rk = Runge_Kutta(f, g, a, y0, z0, h, start);
    for (int k = 1; k <= start; k++)
    {
        s.y[k] = rk.y[k];
        s.z[k] = rk.z[k];
    }

    // значения правых частей в узлах; запоминаем, чтобы не считать их заново на следующих шагах
    Vector fk(N + 1), gk(N + 1);
    for (int k = 0; k <= start; k++)
    {
        fk[k] = f(s.x[k], s.y[k], s.z[k]);
        gk[k] = g(s.x[k], s.y[k], s.z[k]);
    }

    for (int k = 3; k < N; k++)
    {
        s.y[k + 1] = s.y[k] + h / 24 * (55 * fk[k] - 59 * fk[k - 1] + 37 * fk[k - 2] - 9 * fk[k - 3]);
        s.z[k + 1] = s.z[k] + h / 24 * (55 * gk[k] - 59 * gk[k - 1] + 37 * gk[k - 2] - 9 * gk[k - 3]);
        fk[k + 1] = f(s.x[k + 1], s.y[k + 1], s.z[k + 1]);
        gk[k + 1] = g(s.x[k + 1], s.y[k + 1], s.z[k + 1]);
    }
    return s;
}

using Method = Solution (*)(const System_Func &, const System_Func &, double, double, double, double, int);

void Solve_Task_4_1(const std::string &filename)
{
    std::cout << "Задание 4.1. Задача Коши для ОДУ 2-го порядка\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле только шаг h, сама задача из варианта записана в коде выше
    double h;
    if (!(file >> h) || h <= 0)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    double a = PI / 2, b = PI / 2 + 1;
    double y0 = PI / 2, z0 = 1 - PI / 2;
    int N;
    if (!Steps_Count(a, b, h, N))
    {
        std::cout << "Шаг не укладывается на отрезке целое число раз\n";
        return;
    }
    // для Рунге-Ромберга считаем еще раз с шагом 2h, для этого число шагов должно быть четным
    bool rr_ok = N % 2 == 0;

    std::cout << "x^2 y'' - 2x y' + (x^2 + 2) y = 0, y(pi/2) = pi/2, y'(pi/2) = 1 - pi/2\n";
    std::cout << "x in [pi/2, pi/2 + 1] = [" << Format_Number(a) << ", " << Format_Number(b) << "]\n";
    std::cout << "Точное решение y = x cos x + x sin x\n";
    std::cout << "Замена z = y': y' = z, z' = (2xz - (x^2 + 2)y) / x^2\n";
    std::cout << "h = " << h << ", N = " << N << "\n\n";

    const char *names[3] = {"Метод Эйлера", "Метод Рунге-Кутты 4-го порядка", "Метод Адамса 4-го порядка"};
    const int orders[3] = {1, 4, 4}; // порядок точности p
    const Method methods[3] = {Euler, Runge_Kutta, Adams};
    double max_err[3], max_R[3] = {0, 0, 0};

    for (int m = 0; m < 3; m++)
    {
        std::cout << names[m] << ", p = " << orders[m] << "\n";
        Solution s = methods[m](F, G, a, y0, z0, h, N);
        max_err[m] = Print_Solution(s.x, s.y, Exact);

        if (rr_ok)
        {
            Solution s2 = methods[m](F, G, a, y0, z0, 2 * h, N / 2);
            std::cout << "Рунге-Ромберг: решение с шагом 2h = " << 2 * h << ", R = (y_h - y_2h) / (2^p - 1)\n";
            max_R[m] = Print_Runge_Romberg(s.x, s.y, s2.y, orders[m], Exact);
        }
        std::cout << "\n";
    }

    if (!rr_ok)
        std::cout << "Число шагов нечетное, решение с шагом 2h не построить, Рунге-Ромберг не считаем\n\n";

    std::cout << "Максимальная погрешность при h = " << h << ":\n";
    const char *short_names[3] = {"Эйлер", "Рунге-Кутта", "Адамс"};
    for (int m = 0; m < 3; m++)
    {
        std::cout << "  " << short_names[m] << ": по точному решению " << Format_Number(max_err[m], 3);
        if (rr_ok)
            std::cout << ", по Рунге-Ромбергу " << Format_Number(max_R[m], 3);
        std::cout << "\n";
    }
}
