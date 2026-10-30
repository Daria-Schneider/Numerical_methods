// Задание 3.1. Интерполяционные многочлены Лагранжа и Ньютона
// функция y = x*sin(x), два набора узлов а) и б), погрешность считаем в точке X*
#include "interpolation.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

const double PI = std::acos(-1.0);
const int GRID_POINTS = 1000; // на столько частей делим отрезок, когда ищем максимум производной

// интерполируемая функция
double F(double x)
{
    return x * std::sin(x);
}

// k-я производная x*sin(x), нужна только для оценки погрешности (3.9).
// Если продифференцировать несколько раз:
// f' = x*cos(x) + sin(x), f'' = -x*sin(x) + 2cos(x),
// f''' = -x*cos(x) - 3sin(x), f'''' = x*sin(x) - 4cos(x),
// то видно, что f^(k)(x) = x*sin(x + k*pi/2) + k*sin(x + (k-1)*pi/2)
double dF(double x, int k)
{
    return x * std::sin(x + k * PI / 2) + k * std::sin(x + (k - 1) * PI / 2);
}

// w(X) = (X - x_0)(X - x_1)...(X - x_n)
double Omega(const Vector &x, double X)
{
    double w = 1.0;
    for (double xj : x)
        w *= X - xj;
    return w;
}

// w'(x_i) = (x_i - x_0)...(x_i - x_(i-1))(x_i - x_(i+1))...(x_i - x_n),
// то есть то же произведение, но без множителя с j = i
double dOmega(const Vector &x, size_t i)
{
    double w = 1.0;
    for (size_t j = 0; j < x.size(); j++)
        if (j != i)
            w *= x[i] - x[j];
    return w;
}

// Многочлен Лагранжа (3.6): L(X) = sum f_i * w(X) / ((X - x_i) * w'(x_i)).
// w(X) / (X - x_i) - это произведение (X - x_j) по всем j, кроме j = i,
// так и считаем, чтобы не делить на X - x_i (X может совпасть с узлом)
double Lagrange(const Vector &x, const Vector &f, double X)
{
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); i++)
    {
        double term = f[i] / dOmega(x, i);
        for (size_t j = 0; j < x.size(); j++)
            if (j != i)
                term *= X - x[j];
        sum += term;
    }
    return sum;
}

// Таблица разделенных разностей (как табл. 3.1 в методичке).
// diff[0] - значения f_i, diff[1] - разности первого порядка f(x_i, x_(i+1)) и т.д.
// Формула (3.7): f(x_i, ..., x_(i+k)) = (f(x_i, ..., x_(i+k-1)) - f(x_(i+1), ..., x_(i+k))) / (x_i - x_(i+k))
Matrix Divided_Differences(const Vector &x, const Vector &f)
{
    size_t n = x.size();
    Matrix diff(n);
    diff[0] = f;
    for (size_t k = 1; k < n; k++)
        for (size_t i = 0; i + k < n; i++)
            diff[k].push_back((diff[k - 1][i] - diff[k - 1][i + 1]) / (x[i] - x[i + k]));
    return diff;
}

// Многочлен Ньютона (3.8):
// P(X) = f(x_0) + (X - x_0) f(x_0, x_1) + (X - x_0)(X - x_1) f(x_0, x_1, x_2) + ...
// коэффициенты - верхние элементы столбцов таблицы разделенных разностей
double Newton(const Vector &x, const Matrix &diff, double X)
{
    double sum = 0.0;
    double product = 1.0; // (X - x_0)...(X - x_(k-1))
    for (size_t k = 0; k < x.size(); k++)
    {
        sum += diff[k][0] * product;
        product *= X - x[k];
    }
    return sum;
}

// Узлы в файле записаны в долях pi, как в условии: "5/12" значит 5*pi/12, "0" значит 0
static bool Parse_Pi(const std::string &s, double &value)
{
    std::istringstream ss(s);
    double num, den = 1.0;
    char slash;
    if (!(ss >> num))
        return false;
    if (ss >> slash && (slash != '/' || !(ss >> den) || den == 0))
        return false;
    value = num / den * PI;
    return true;
}

// для вывода: "5/12" -> "5pi/12", "1/6" -> "pi/6", "0" -> "0"
static std::string Pi_Text(const std::string &s)
{
    size_t slash = s.find('/');
    std::string num = s.substr(0, slash);
    std::string rest = (slash == std::string::npos) ? "" : s.substr(slash);
    if (num == "0")
        return "0";
    if (num == "1")
        num = "";
    return num + "pi" + rest;
}

// строим оба многочлена для одного набора узлов, возвращаем погрешность в X*
static double Solve_For_Nodes(const Vector &x, double X_star)
{
    size_t n = x.size();
    Vector f(n);
    for (size_t i = 0; i < n; i++)
        f[i] = F(x[i]);

    // Лагранж: таблица как в примере 3.1 из методички
    std::cout << "Многочлен Лагранжа степени " << n - 1 << "\n";
    std::cout << std::setw(3) << "i" << std::setw(12) << "x_i" << std::setw(12) << "f_i" << std::setw(12)
              << "w'(x_i)" << std::setw(14) << "f_i/w'(x_i)" << std::setw(12) << "X* - x_i" << "\n";

    std::string lagrange_text;
    for (size_t i = 0; i < n; i++)
    {
        double dw = dOmega(x, i);
        std::cout << std::setw(3) << i << std::setw(12) << Format_Number(x[i]) << std::setw(12)
                  << Format_Number(f[i]) << std::setw(12) << Format_Number(dw) << std::setw(14)
                  << Format_Number(f[i] / dw) << std::setw(12) << Format_Number(X_star - x[i]) << "\n";

        std::string factors;
        for (size_t j = 0; j < n; j++)
            if (j != i)
                factors += Factor_Text(x[j]);
        Add_Term(lagrange_text, f[i] / dw, factors);
    }
    std::cout << "L(x) = " << (lagrange_text.empty() ? "0" : lagrange_text) << "\n\n";

    // Ньютон: таблица разделенных разностей, в строке i стоят f_i, f(x_i, x_(i+1)), ...
    Matrix diff = Divided_Differences(x, f);
    std::cout << "Многочлен Ньютона степени " << n - 1 << ", таблица разделенных разностей:\n";
    std::cout << std::setw(3) << "i" << std::setw(12) << "x_i" << std::setw(15) << "f_i";
    for (size_t k = 1; k < n; k++)
        std::cout << std::setw(15) << (k == 1 ? "f(x_i,x_i+1)" : "f(x_i..x_i+" + std::to_string(k) + ")");
    std::cout << "\n";
    for (size_t i = 0; i < n; i++)
    {
        std::cout << std::setw(3) << i << std::setw(12) << Format_Number(x[i]);
        for (size_t k = 0; i + k < n; k++)
            std::cout << std::setw(15) << Format_Number(diff[k][i]);
        std::cout << "\n";
    }

    std::string newton_text, factors;
    for (size_t k = 0; k < n; k++)
    {
        Add_Term(newton_text, diff[k][0], factors);
        factors += Factor_Text(x[k]);
    }
    std::cout << "P(x) = " << (newton_text.empty() ? "0" : newton_text) << "\n\n";

    // значения в точке X*
    double L = Lagrange(x, f, X_star);
    double P = Newton(x, diff, X_star);
    double y = F(X_star);
    double error = std::abs(y - L);
    std::cout << "L(X*) = " << Format_Number(L) << ", P(X*) = " << Format_Number(P) << "\n";
    std::cout << "Точное значение y(X*) = " << Format_Number(y) << "\n";
    std::cout << "Погрешность интерполяции |y(X*) - L(X*)| = " << Format_Number(error, 4) << "\n";

    // Проверка по оценке (3.9): |f(X*) - L(X*)| <= M / (n+1)! * |w(X*)|,
    // где M = max|f^(n+1)| на [x_0, x_n] (n - степень многочлена, у нас узлов n + 1)
    double a = *std::min_element(x.begin(), x.end());
    double b = *std::max_element(x.begin(), x.end());
    int order = static_cast<int>(n);
    double M = 0.0;
    for (int i = 0; i <= GRID_POINTS; i++)
        M = std::max(M, std::abs(dF(a + (b - a) * i / GRID_POINTS, order)));

    double factorial = 1.0;
    for (int k = 2; k <= order; k++)
        factorial *= k;
    double w = std::abs(Omega(x, X_star));
    double estimate = M / factorial * w;

    std::cout << "Оценка (3.9): M = max|f^(" << order << ")(x)| = " << Format_Number(M) << ", |w(X*)| = "
              << Format_Number(w) << ", M / " << order << "! * |w(X*)| = " << Format_Number(estimate, 4) << "\n";
    std::cout << "Погрешность не больше оценки: " << (error <= estimate ? "верно" : "неверно") << "\n\n";

    return error;
}

void Solve_Task_3_1(const std::string &filename)
{
    std::cout << "Задание 3.1. Интерполяционные многочлены Лагранжа и Ньютона\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле три строки: узлы для а), узлы для б) и точка X*, все в долях pi
    std::vector<std::string> lines[3];
    for (auto &line : lines)
    {
        if (!Read_Tokens(file, line))
        {
            std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
            return;
        }
    }

    Vector nodes[2];
    double X_star;
    bool ok = Parse_Pi(lines[2][0], X_star);
    for (int v = 0; v < 2; v++)
    {
        for (const std::string &s : lines[v])
        {
            double value;
            ok = ok && Parse_Pi(s, value);
            nodes[v].push_back(value);
        }
        ok = ok && nodes[v].size() >= 2;
        // узлы не должны совпадать, иначе в w'(x_i) будет ноль
        for (size_t i = 0; ok && i < nodes[v].size(); i++)
            for (size_t j = i + 1; j < nodes[v].size(); j++)
                if (std::abs(nodes[v][i] - nodes[v][j]) < 1e-12)
                    ok = false;
    }
    if (!ok)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    std::cout << "Функция y = x*sin(x), X* = " << Pi_Text(lines[2][0]) << " = " << Format_Number(X_star) << "\n\n";

    const char *names[2] = {"а)", "б)"};
    double errors[2];
    for (int v = 0; v < 2; v++)
    {
        std::cout << names[v] << " X_i = ";
        for (size_t i = 0; i < lines[v].size(); i++)
            std::cout << (i > 0 ? ", " : "") << Pi_Text(lines[v][i]);
        std::cout << "\n\n";

        errors[v] = Solve_For_Nodes(nodes[v], X_star);
    }

    std::cout << "Погрешность в точке X*: а) " << Format_Number(errors[0], 4) << ", б) "
              << Format_Number(errors[1], 4) << "\n";
}
