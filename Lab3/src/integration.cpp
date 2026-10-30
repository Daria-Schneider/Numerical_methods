// Задание 3.5. Численное интегрирование: методы прямоугольников, трапеций и Симпсона,
// оценка погрешности методом Рунге-Ромберга
#include "integration.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// подынтегральная функция
double Y(double x)
{
    return std::sqrt(9 + x * x);
}

// Точное значение, только для проверки. Первообразная (табличный интеграл):
// integral sqrt(9 + x^2) dx = x/2 * sqrt(9 + x^2) + 9/2 * ln(x + sqrt(9 + x^2))
double Exact_Integral(double a, double b)
{
    auto G = [](double x) { return x / 2 * std::sqrt(9 + x * x) + 4.5 * std::log(x + std::sqrt(9 + x * x)); };
    return G(b) - G(a);
}

// Все три функции ниже возвращают интеграл от a до каждого узла x_i = a + i*h,
// как в таблицах методички. Последний элемент - интеграл по всему отрезку.

// Метод прямоугольников (3.23): на каждом отрезке берем значение функции в середине
Vector Rectangles(double a, double h, int N)
{
    Vector res(N + 1, 0.0);
    for (int i = 1; i <= N; i++)
        res[i] = res[i - 1] + h * Y(a + (i - 0.5) * h);
    return res;
}

// Метод трапеций (3.25): на каждом отрезке площадь трапеции h * (y_(i-1) + y_i) / 2
Vector Trapezoids(double a, double h, int N)
{
    Vector res(N + 1, 0.0);
    for (int i = 1; i <= N; i++)
        res[i] = res[i - 1] + h * (Y(a + (i - 1) * h) + Y(a + i * h)) / 2;
    return res;
}

// Метод Симпсона (3.28): через каждые три точки x_(2m-2), x_(2m-1), x_(2m) проводим параболу,
// ее интеграл h/3 * (y_(2m-2) + 4y_(2m-1) + y_(2m)). Поэтому N должно быть четным,
// а значения получаются только в узлах с четными номерами: res[m] - интеграл до x_(2m)
Vector Simpson(double a, double h, int N)
{
    Vector res(N / 2 + 1, 0.0);
    for (int m = 1; m <= N / 2; m++)
    {
        double y0 = Y(a + (2 * m - 2) * h);
        double y1 = Y(a + (2 * m - 1) * h);
        double y2 = Y(a + 2 * m * h);
        res[m] = res[m - 1] + h / 3 * (y0 + 4 * y1 + y2);
    }
    return res;
}

// Метод Рунге-Ромберга (3.30): если F_h посчитан с шагом h, а F_kh с шагом kh,
// и у метода порядок точности p, то F = F_h + (F_h - F_kh) / (k^p - 1) + O(h^(p+1)).
// Добавка (F_h - F_kh) / (k^p - 1) - это и есть оценка погрешности F_h
double Runge_Romberg_Error(double F_h, double F_kh, double k, int p)
{
    return (F_h - F_kh) / (std::pow(k, p) - 1);
}

// таблица для одного шага, возвращает значения по всему отрезку для трех методов
static void Print_Table(double a, double h, int N, bool simpson_ok, double &rect, double &trap, double &simp)
{
    Vector r = Rectangles(a, h, N);
    Vector t = Trapezoids(a, h, N);
    Vector s;
    if (simpson_ok)
        s = Simpson(a, h, N);

    std::cout << "Шаг h = " << h << ", отрезков N = " << N << "\n";
    std::cout << "в столбцах методов - интеграл от X0 до x_i\n";
    std::cout << "   i       x_i         y_i  прямоугольники    трапеции     Симпсон\n";
    for (int i = 0; i <= N; i++)
    {
        double x = a + i * h;
        std::cout << std::setw(4) << i << std::setw(10) << Format_Number(x) << std::setw(12) << Format_Number(Y(x))
                  << std::setw(16) << Format_Number(r[i], 7) << std::setw(12) << Format_Number(t[i], 7);
        if (simpson_ok && i % 2 == 0)
            std::cout << std::setw(12) << Format_Number(s[i / 2], 7);
        else
            std::cout << std::setw(12) << "-";
        std::cout << "\n";
    }
    std::cout << "\n";

    rect = r[N];
    trap = t[N];
    simp = simpson_ok ? s[N / 2] : 0.0;
}

void Solve_Task_3_5(const std::string &filename)
{
    std::cout << "Задание 3.5. Численное интегрирование\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле: пределы X0, Xk и шаги h1, h2 (h1 > h2)
    double a, b, h1, h2;
    if (!(file >> a >> b >> h1 >> h2) || a >= b || h2 <= 0 || h1 <= h2)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    // шаг должен укладываться на отрезке целое число раз
    int N1 = static_cast<int>(std::round((b - a) / h1));
    int N2 = static_cast<int>(std::round((b - a) / h2));
    if (std::abs(N1 * h1 - (b - a)) > 1e-9 || std::abs(N2 * h2 - (b - a)) > 1e-9)
    {
        std::cout << "Шаг не укладывается на отрезке целое число раз\n";
        return;
    }
    bool simpson_ok = N1 % 2 == 0 && N2 % 2 == 0;

    std::cout << "y = sqrt(9 + x^2), X0 = " << a << ", Xk = " << b << ", h1 = " << h1 << ", h2 = " << h2 << "\n\n";

    double F1[3], F2[3]; // результаты с шагом h1 и h2: прямоугольники, трапеции, Симпсон
    Print_Table(a, h1, N1, simpson_ok, F1[0], F1[1], F1[2]);
    Print_Table(a, h2, N2, simpson_ok, F2[0], F2[1], F2[2]);
    if (!simpson_ok)
        std::cout << "Для метода Симпсона число отрезков должно быть четным, его не считаем\n\n";

    double exact = Exact_Integral(a, b);
    std::cout << "Точное значение F = " << Format_Number(exact, 8) << "\n\n";

    // Рунге-Ромберг: h = h2 (мелкий шаг), kh = h1, значит k = h1 / h2
    double k = h1 / h2;
    std::cout << "Оценка погрешности методом Рунге-Ромберга, k = h1 / h2 = " << Format_Number(k) << ":\n";
    std::cout << "R = (F(h2) - F(h1)) / (k^p - 1), уточненное значение F(h2) + R\n\n";

    const char *names[3] = {"Метод прямоугольников", "Метод трапеций", "Метод Симпсона"};
    const int orders[3] = {2, 2, 4}; // порядок точности p
    for (int m = 0; m < 3; m++)
    {
        if (m == 2 && !simpson_ok)
            break;

        double R = Runge_Romberg_Error(F2[m], F1[m], k, orders[m]);
        double refined = F2[m] + R;
        std::cout << names[m] << ", p = " << orders[m] << ":\n";
        std::cout << "  F(h1) = " << Format_Number(F1[m], 8) << ", F(h2) = " << Format_Number(F2[m], 8) << "\n";
        std::cout << "  R = " << Format_Number(R, 4) << ", уточненное значение " << Format_Number(refined, 8) << "\n";
        std::cout << "  отличие от точного значения: F(h1): " << Format_Number(std::abs(exact - F1[m]), 3)
                  << ", F(h2): " << Format_Number(std::abs(exact - F2[m]), 3)
                  << ", уточненное: " << Format_Number(std::abs(exact - refined), 3) << "\n\n";
    }
}
