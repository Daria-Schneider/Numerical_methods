// Задание 3.4. Численное дифференцирование таблично заданной функции
#include "differentiation.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// (3.18): функцию на [x_i, x_(i+1)] заменяем прямой, ее наклон и есть производная.
// Первый порядок точности
double First_Derivative_1(const Vector &x, const Vector &y, size_t i)
{
    return (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

// (3.20): функцию заменяем многочленом 2-й степени по точкам x_i, x_(i+1), x_(i+2)
// (это многочлен Ньютона из (3.19)) и дифференцируем его:
// y'(X) = (y_(i+1) - y_i) / (x_(i+1) - x_i) + [вторая разделенная разность] * (2X - x_i - x_(i+1)).
// На равномерной сетке второй порядок точности
double First_Derivative_2(const Vector &x, const Vector &y, size_t i, double X)
{
    double left = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    double right = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
    double second_diff = (right - left) / (x[i + 2] - x[i]);
    return left + second_diff * (2 * X - x[i] - x[i + 1]);
}

// (3.21): вторая производная того же многочлена 2-й степени,
// y'' = 2 * [вторая разделенная разность], она постоянная на отрезке
double Second_Derivative(const Vector &x, const Vector &y, size_t i)
{
    double left = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    double right = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
    return 2 * (right - left) / (x[i + 2] - x[i]);
}

void Solve_Task_3_4(const std::string &filename)
{
    std::cout << "Задание 3.4. Численное дифференцирование\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле: точка X*, строка с x_i, строка с y_i
    Vector star, x, y;
    bool ok = Read_Numbers(file, star) && Read_Numbers(file, x) && Read_Numbers(file, y);
    ok = ok && x.size() == y.size() && x.size() >= 3;
    for (size_t i = 1; ok && i < x.size(); i++)
        if (x[i] <= x[i - 1])
            ok = false;
    if (!ok)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }
    double X = star[0];
    size_t n = x.size() - 1;

    std::cout << std::setw(5) << "i";
    for (size_t i = 0; i <= n; i++)
        std::cout << std::setw(10) << i;
    std::cout << "\n" << std::setw(5) << "x_i";
    for (size_t i = 0; i <= n; i++)
        std::cout << std::setw(10) << Format_Number(x[i]);
    std::cout << "\n" << std::setw(5) << "y_i";
    for (size_t i = 0; i <= n; i++)
        std::cout << std::setw(10) << Format_Number(y[i]);
    std::cout << "\n\n";

    // ищем отрезок [x_k, x_(k+1)], где лежит X*.
    // Если X* совпадает с узлом, берется отрезок, у которого X* - правый конец (как в примере из методички)
    size_t k = n;
    for (size_t i = 0; i < n; i++)
    {
        if (x[i] <= X && X <= x[i + 1])
        {
            k = i;
            break;
        }
    }
    if (k == n)
    {
        std::cout << "X* = " << X << " не лежит на отрезке [" << x[0] << ", " << x[n] << "]\n";
        return;
    }
    bool is_node = std::abs(X - x[k + 1]) < 1e-12; // X* совпадает с узлом x_(k+1)

    std::cout << "X* = " << X << ", лежит на отрезке [x_" << k << ", x_" << k + 1 << "] = [" << x[k] << ", "
              << x[k + 1] << "]\n\n";

    std::cout << "Первая производная с первым порядком точности (3.18):\n";
    double d_left = First_Derivative_1(x, y, k);
    std::string side = "";
    if (is_node)
        side = "левосторонняя, ";
    else if (std::abs(X - x[k]) < 1e-12) // X* = x_0, левее точек нет
        side = "правосторонняя, ";
    std::cout << "  " << side << "отрезок [" << x[k] << ", " << x[k + 1] << "]: y'(X*) = (y_" << k + 1 << " - y_" << k
              << ") / (x_" << k + 1 << " - x_" << k << ") = " << Format_Number(d_left) << "\n";

    // если X* - узел и справа есть еще точка, считаем и правостороннюю производную
    bool has_right = is_node && k + 2 <= n;
    double d_right = 0.0;
    if (has_right)
    {
        d_right = First_Derivative_1(x, y, k + 1);
        std::cout << "  правосторонняя, отрезок [" << x[k + 1] << ", " << x[k + 2] << "]: y'(X*) = (y_" << k + 2
                  << " - y_" << k + 1 << ") / (x_" << k + 2 << " - x_" << k + 1 << ") = " << Format_Number(d_right)
                  << "\n";
    }
    std::cout << "\n";

    // для (3.20) и (3.21) нужны три точки x_i, x_(i+1), x_(i+2), у правого края сдвигаемся влево
    size_t i = std::min(k, n - 2);
    std::cout << "Многочлен 2-й степени по точкам x_" << i << ", x_" << i + 1 << ", x_" << i + 2 << " = " << x[i]
              << ", " << x[i + 1] << ", " << x[i + 2] << "\n";

    double d1 = First_Derivative_2(x, y, i, X);
    std::cout << "Первая производная со вторым порядком точности (3.20): y'(X*) = " << Format_Number(d1) << "\n";
    if (has_right)
        std::cout << "Полусумма левой и правой производных: " << Format_Number((d_left + d_right) / 2)
                  << " (на равномерной сетке должна совпасть с (3.20))\n";

    double d2 = Second_Derivative(x, y, i);
    std::cout << "Вторая производная (3.21): y''(X*) = " << Format_Number(d2) << "\n";
}
