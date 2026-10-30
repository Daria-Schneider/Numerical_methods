// Задание 2.1. Метод простой итерации и метод Ньютона
// уравнение x^6 - 5x^3 - 2 = 0, ищем положительный корень
#include "nonlinear_equation.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

const int MAX_ITERATIONS = 1000;
const int GRID_POINTS = 1000; // на столько частей делим отрезок, когда проверяем условия

// f(x) - левая часть уравнения
double F(double x)
{
    return std::pow(x, 6) - 5 * std::pow(x, 3) - 2;
}

// f'(x)
double dF(double x)
{
    return 6 * std::pow(x, 5) - 15 * x * x;
}

// f''(x), нужна только чтобы правильно выбрать x0 для метода Ньютона
double d2F(double x)
{
    return 30 * std::pow(x, 4) - 30 * x;
}

// Для простой итерации уравнение надо записать в виде x = phi(x):
// x^6 = 5x^3 + 2  ->  x = (5x^3 + 2)^(1/6)
double Phi(double x)
{
    return std::pow(5 * std::pow(x, 3) + 2, 1.0 / 6);
}

// phi'(x) = 1/6 * (5x^3 + 2)^(-5/6) * 15x^2 = 2.5x^2 * (5x^3 + 2)^(-5/6)
double dPhi(double x)
{
    return 2.5 * x * x * std::pow(5 * std::pow(x, 3) + 2, -5.0 / 6);
}

// Точный корень, нужен только для проверки.
// Замена t = x^3 дает квадратное уравнение t^2 - 5t - 2 = 0,
// его положительный корень t = (5 + sqrt(33)) / 2, значит x = t^(1/3)
double Exact_Root()
{
    return std::cbrt((5 + std::sqrt(33.0)) / 2);
}

// q = max|phi'(x)| на [a, b], просто перебираем точки отрезка
double Find_q(double a, double b)
{
    double q = 0.0;
    for (int i = 0; i <= GRID_POINTS; i++)
    {
        double x = a + (b - a) * i / GRID_POINTS;
        q = std::max(q, std::abs(dPhi(x)));
    }
    return q;
}

// наименьшее и наибольшее значение phi(x) на [a, b]
// (по теореме 2.3 phi(x) не должна выходить за отрезок)
void Phi_Range(double a, double b, double &phi_min, double &phi_max)
{
    phi_min = Phi(a);
    phi_max = Phi(a);
    for (int i = 0; i <= GRID_POINTS; i++)
    {
        double x = a + (b - a) * i / GRID_POINTS;
        phi_min = std::min(phi_min, Phi(x));
        phi_max = std::max(phi_max, Phi(x));
    }
}

// true, если функция g на [a, b] не меняет знак (проверяем f' и f'' для Ньютона)
bool Keeps_Sign(double (*g)(double), double a, double b)
{
    double first = g(a);
    for (int i = 0; i <= GRID_POINTS; i++)
    {
        double x = a + (b - a) * i / GRID_POINTS;
        if (g(x) * first <= 0)
            return false;
    }
    return true;
}

// Метод простой итерации: x(k+1) = phi(x(k)).
// Останавливаемся по оценке (2.7): |x* - x(k+1)| <= q / (1 - q) * |x(k+1) - x(k)|,
// coef - это и есть q / (1 - q)
Equation_Result Simple_Iteration(double x0, double eps, double coef)
{
    Equation_Result res;
    double x = x0;
    res.history.push_back(x);

    for (int k = 1; k <= MAX_ITERATIONS; k++)
    {
        double x_new = Phi(x);
        double err = coef * std::abs(x_new - x); // оценка погрешности
        x = x_new;

        res.iterations = k;
        res.history.push_back(x);
        res.errors.push_back(err);

        if (!std::isfinite(err)) // расходится
            break;
        if (err <= eps)
        {
            res.converged = true;
            break;
        }
    }
    res.x = x;
    return res;
}

// Метод Ньютона (касательных): x(k+1) = x(k) - f(x(k)) / f'(x(k)).
// Останавливаемся, когда |x(k+1) - x(k)| <= eps
Equation_Result Newton_Method(double x0, double eps)
{
    Equation_Result res;
    double x = x0;
    res.history.push_back(x);

    for (int k = 1; k <= MAX_ITERATIONS; k++)
    {
        double d = dF(x);
        if (d == 0) // касательная горизонтальная, на 0 делить нельзя
            break;

        double x_new = x - F(x) / d;
        double err = std::abs(x_new - x);
        x = x_new;

        res.iterations = k;
        res.history.push_back(x);
        res.errors.push_back(err);

        if (!std::isfinite(err))
            break;
        if (err <= eps)
        {
            res.converged = true;
            break;
        }
    }
    res.x = x;
    return res;
}

// таблица итераций: x(k), f(x(k)), оценка погрешности и настоящая погрешность |x(k) - x*|
static void Print_Iterations(const Equation_Result &res, double x_star)
{
    std::cout << std::setw(4) << "k" << std::setw(13) << "x(k)" << std::setw(13) << "f(x(k))" << std::setw(13)
              << "eps(k)" << std::setw(15) << "|x(k) - x*|" << "\n";

    for (size_t k = 0; k < res.history.size(); k++)
    {
        double x = res.history[k];
        std::cout << std::setw(4) << k << std::setw(13) << Format_Number(x) << std::setw(13) << Format_Number(F(x));
        if (k == 0)
            std::cout << std::setw(13) << "-";
        else
            std::cout << std::setw(13) << Format_Number(res.errors[k - 1], 3);
        std::cout << std::setw(15) << Format_Number(std::abs(x - x_star), 3) << "\n";
    }
    std::cout << "\n";
}

void Solve_Task_2_1(const std::string &filename)
{
    std::cout << "Задание 2.1. Метод простой итерации и метод Ньютона\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле: точность eps и отрезок [a, b], на котором по графику лежит корень
    double eps, a, b;
    if (!(file >> eps >> a >> b) || eps <= 0 || a >= b)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    std::cout << "Уравнение: x^6 - 5x^3 - 2 = 0\n";
    std::cout << "Точность eps = " << eps << "\n";
    std::cout << "Отрезок с положительным корнем (по графику): [" << a << ", " << b << "]\n";
    std::cout << "f(" << a << ") = " << Format_Number(F(a)) << ", f(" << b << ") = " << Format_Number(F(b)) << "\n";

    if (F(a) * F(b) >= 0)
    {
        std::cout << "На концах отрезка f(x) одного знака, нужно взять другой отрезок\n";
        return;
    }
    std::cout << "f(a) * f(b) < 0, значит на отрезке есть корень\n\n";

    // Метод Ньютона, проверяем условия теоремы 2.2
    std::cout << "Метод Ньютона\n";
    std::cout << "f'(x) = 6x^5 - 15x^2, f''(x) = 30x^4 - 30x\n";
    bool signs_ok = Keeps_Sign(dF, a, b) && Keeps_Sign(d2F, a, b);
    std::cout << "f'(x) и f''(x) на отрезке знак " << (signs_ok ? "не меняют" : "меняют") << "\n";

    // x0 берем на том конце отрезка, где f(x0) * f''(x0) > 0 (условие 2.3),
    // тогда метод Ньютона сходится к корню монотонно
    double x0 = b;
    if (F(a) * d2F(a) > 0)
        x0 = a;
    std::cout << "x0 = " << x0 << ", f(x0) * f''(x0) = " << Format_Number(F(x0) * d2F(x0));
    if (F(x0) * d2F(x0) > 0)
        std::cout << " > 0\n\n";
    else
        std::cout << " <= 0, сходимость не гарантирована\n\n";

    // Метод простой итерации, проверяем условия теоремы 2.3
    std::cout << "Метод простой итерации\n";
    std::cout << "phi(x) = (5x^3 + 2)^(1/6)\n";

    double phi_min, phi_max;
    Phi_Range(a, b, phi_min, phi_max);
    bool phi_inside = phi_min >= a && phi_max <= b;
    std::cout << "phi(x) на отрезке принимает значения от " << Format_Number(phi_min) << " до "
              << Format_Number(phi_max) << (phi_inside ? ", за отрезок не выходит\n" : ", выходит за отрезок\n");

    double q = Find_q(a, b);
    std::cout << "q = max|phi'(x)| = " << Format_Number(q) << "\n";

    // если условия выполнены, погрешность оцениваем как q / (1 - q) * |x(k) - x(k-1)|,
    // иначе просто |x(k) - x(k-1)|
    double coef = 1.0;
    bool simple_ok = phi_inside && q < 1;
    if (simple_ok)
    {
        coef = q / (1 - q);
        std::cout << "q < 1, условия теоремы 2.3 выполнены, метод сходится\n";
    }
    else
    {
        std::cout << "Условия теоремы 2.3 не выполнены, сходимость не гарантирована\n";
    }
    std::cout << "x0 = " << x0 << " (как у метода Ньютона, чтобы сравнивать методы)\n\n";

    double x_star = Exact_Root();

    Equation_Result simple = Simple_Iteration(x0, eps, coef);
    std::cout << "Метод простой итерации, eps(k) = " << (simple_ok ? "q / (1 - q) * " : "")
              << "|x(k) - x(k-1)|:\n";
    Print_Iterations(simple, x_star);

    Equation_Result newton = Newton_Method(x0, eps);
    std::cout << "Метод Ньютона, eps(k) = |x(k) - x(k-1)|:\n";
    Print_Iterations(newton, x_star);

    std::cout << "Корень методом простой итерации: x = " << Format_Number(simple.x)
              << ", итераций: " << simple.iterations << (simple.converged ? "" : " (не сошелся)") << "\n";
    std::cout << "Корень методом Ньютона: x = " << Format_Number(newton.x) << ", итераций: " << newton.iterations
              << (newton.converged ? "" : " (не сошелся)") << "\n";
    std::cout << "Точный корень (через замену t = x^3): x* = " << Format_Number(x_star, 10) << "\n";

    // как число итераций зависит от точности
    std::cout << "\nКоличество итераций при разной точности:\n";
    std::cout << std::setw(10) << "eps" << "   простая итерация   Ньютон\n";
    for (double e : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6})
    {
        Equation_Result s = Simple_Iteration(x0, e, coef);
        Equation_Result n = Newton_Method(x0, e);
        std::cout << std::setw(10) << e << std::setw(19) << s.iterations << std::setw(9) << n.iterations << "\n";
    }
    std::cout << "\n";

    if (simple.converged && newton.converged)
    {
        if (newton.iterations < simple.iterations)
            std::cout << "Метод Ньютона сошелся быстрее: итераций " << newton.iterations << " против "
                      << simple.iterations << "\n";
        else if (newton.iterations > simple.iterations)
            std::cout << "Метод простой итерации сошелся быстрее: итераций " << simple.iterations << " против "
                      << newton.iterations << "\n";
        else
            std::cout << "Оба метода сошлись за одинаковое число итераций\n";
    }
}
