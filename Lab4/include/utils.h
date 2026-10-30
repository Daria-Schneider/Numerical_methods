#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <string>
#include <vector>

using Vector = std::vector<double>;
using Func = std::function<double(double)>;

// число для вывода в таблицы, без лишних нулей
std::string Format_Number(double value, int precision = 6);

// число шагов N = (b - a) / h, false если шаг не укладывается на отрезке целое число раз
bool Steps_Count(double a, double b, double h, int &N);

// оценка погрешности решения с шагом h по решениям с шагами h и 2h
double Runge_Romberg_Error(double y_h, double y_2h, int p);

// таблица k, x_k, y_k, точное значение, погрешность. Возвращает максимальную погрешность
double Print_Solution(const Vector &x, const Vector &y, const Func &exact);

// таблица метода Рунге-Ромберга в общих узлах сеток h и 2h. Возвращает максимальное |R|
double Print_Runge_Romberg(const Vector &x, const Vector &y_h, const Vector &y_2h, int p, const Func &exact);

#endif // UTILS_H
