#include "utils.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

std::string Format_Number(double value, int precision)
{
    // погрешности округления вида 1e-17 печатаем как 0
    if (std::abs(value) < 1e-10)
        value = 0.0;

    std::ostringstream oss;
    oss << std::setprecision(precision) << value;
    return oss.str();
}

bool Steps_Count(double a, double b, double h, int &N)
{
    N = static_cast<int>(std::round((b - a) / h));
    return N > 0 && std::abs(N * h - (b - a)) < 1e-9;
}

// Метод Рунге-Ромберга (4.11). y_h и y_2h - решения в одной и той же точке, посчитанные одним методом
// с шагами h и 2h, p - порядок точности метода. Погрешность метода примерно C*h^p, то есть
//   y - y_h ~ C h^p,   y - y_2h ~ C (2h)^p = 2^p C h^p.
// Вычитаем одно из другого: y_h - y_2h ~ (2^p - 1) C h^p, отсюда погрешность
//   y - y_h ~ (y_h - y_2h) / (2^p - 1).
// Точное решение для этого не нужно, поэтому так оценивают погрешность на практике
double Runge_Romberg_Error(double y_h, double y_2h, int p)
{
    return (y_h - y_2h) / (std::pow(2.0, p) - 1);
}

double Print_Solution(const Vector &x, const Vector &y, const Func &exact)
{
    double max_err = 0;
    std::cout << "   k       x_k           y_k      y точное   погрешность\n";
    for (size_t k = 0; k < x.size(); k++)
    {
        double err = std::abs(exact(x[k]) - y[k]);
        max_err = std::max(max_err, err);
        std::cout << std::setw(4) << k << std::setw(10) << Format_Number(x[k]) << std::setw(14)
                  << Format_Number(y[k], 8) << std::setw(14) << Format_Number(exact(x[k]), 8) << std::setw(14)
                  << Format_Number(err, 3) << "\n";
    }
    return max_err;
}

double Print_Runge_Romberg(const Vector &x, const Vector &y_h, const Vector &y_2h, int p, const Func &exact)
{
    double max_R = 0;
    std::cout << "       x           y_h          y_2h             R       y - y_h\n";
    for (size_t j = 0; j < y_2h.size(); j++)
    {
        size_t k = 2 * j; // j-й узел сетки 2h совпадает с (2j)-м узлом сетки h
        double R = Runge_Romberg_Error(y_h[k], y_2h[j], p);
        max_R = std::max(max_R, std::abs(R));
        std::cout << std::setw(8) << Format_Number(x[k]) << std::setw(14) << Format_Number(y_h[k], 8)
                  << std::setw(14) << Format_Number(y_2h[j], 8) << std::setw(14) << Format_Number(R, 3)
                  << std::setw(14) << Format_Number(exact(x[k]) - y_h[k], 3) << "\n";
    }
    return max_R;
}
