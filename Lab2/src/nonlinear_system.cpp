// Задание 2.2. Метод простой итерации и метод Ньютона для системы
//   f1 = a*x1^2 - x2 + x2^2 - a = 0
//   f2 = x1 - sqrt(x2 + a) + 1 = 0
// ищем решение, в котором x1 > 0 и x2 > 0
#include "nonlinear_system.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

const int MAX_ITERATIONS = 1000;
const int GRID_POINTS = 100; // столько шагов по каждой оси, когда ищем q в области G

// значения f1 и f2 в точке x = (x1, x2)
Vector F_System(const Vector &x, double a)
{
    double f1 = a * x[0] * x[0] - x[1] + x[1] * x[1] - a;
    double f2 = x[0] - std::sqrt(x[1] + a) + 1;
    return {f1, f2};
}

// матрица Якоби J(x), в строке i производные f_i по x1 и по x2:
// df1/dx1 = 2a*x1    df1/dx2 = 2x2 - 1
// df2/dx1 = 1        df2/dx2 = -1 / (2 * sqrt(x2 + a))
Matrix Jacobian(const Vector &x, double a)
{
    return {{2 * a * x[0], 2 * x[1] - 1},
            {1.0, -1.0 / (2 * std::sqrt(x[1] + a))}};
}

double Det2(const Matrix &M)
{
    return M[0][0] * M[1][1] - M[0][1] * M[1][0];
}

// обратная матрица 2x2: элементы главной диагонали меняем местами,
// у побочной меняем знак и все делим на определитель
Matrix Inverse2(const Matrix &M)
{
    double det = Det2(M);
    return {{M[1][1] / det, -M[0][1] / det},
            {-M[1][0] / det, M[0][0] / det}};
}

// ||x - y|| = max|x_i - y_i|
double Diff_Norm(const Vector &x, const Vector &y)
{
    return std::max(std::abs(x[0] - y[0]), std::abs(x[1] - y[1]));
}

// норма матрицы - максимальная по строкам сумма модулей, как в формуле (2.27)
double Matrix_Norm(const Matrix &M)
{
    double norm = 0.0;
    for (const auto &row : M)
        norm = std::max(norm, std::abs(row[0]) + std::abs(row[1]));
    return norm;
}

bool In_Region(const Vector &x, const Region &G)
{
    return x[0] >= G.x1_min && x[0] <= G.x1_max && x[1] >= G.x2_min && x[1] <= G.x2_max;
}

// Для простой итерации систему надо записать в виде x = phi(x).
// Если просто выразить x2 из одного из уравнений, например x2 = (x1 + 1)^2 - a,
// то |dphi2/dx1| около корня получается от 3 до 6, это больше 1 и условие (2.25) не выполняется.
// Поэтому, как в замечании из методички (x = x - lambda * f(x)), берем
//   phi(x) = x - L * f(x),  где L = J^(-1)(x0) - матрица, обратная к матрице Якоби в x0.
// Корни у x = phi(x) те же, что у f(x) = 0, потому что L невырожденная.
Vector Phi_System(const Vector &x, const Matrix &L, double a)
{
    Vector f = F_System(x, a);
    return {x[0] - (L[0][0] * f[0] + L[0][1] * f[1]),
            x[1] - (L[1][0] * f[0] + L[1][1] * f[1])};
}

// phi'(x) = E - L * J(x). Около x0 произведение L * J(x) почти единичная матрица,
// поэтому phi'(x) маленькая и q < 1
Matrix dPhi_System(const Vector &x, const Matrix &L, double a)
{
    Matrix J = Jacobian(x, a);
    Matrix D(2, Vector(2, 0.0));
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            D[i][j] = (i == j ? 1.0 : 0.0) - (L[i][0] * J[0][j] + L[i][1] * J[1][j]);
    return D;
}

// q = max ||phi'(x)|| по области G, перебираем точки сетки
double Find_q_System(const Matrix &L, double a, const Region &G)
{
    double q = 0.0;
    for (int i = 0; i <= GRID_POINTS; i++)
    {
        for (int j = 0; j <= GRID_POINTS; j++)
        {
            Vector x = {G.x1_min + (G.x1_max - G.x1_min) * i / GRID_POINTS,
                        G.x2_min + (G.x2_max - G.x2_min) * j / GRID_POINTS};
            q = std::max(q, Matrix_Norm(dPhi_System(x, L, a)));
        }
    }
    return q;
}

// Метод простой итерации: x(k+1) = phi(x(k)).
// Останавливаемся по оценке (2.26): q / (1 - q) * ||x(k+1) - x(k)|| <= eps,
// coef - это q / (1 - q)
System_Result Simple_Iteration_System(const Vector &x0, const Matrix &L, double a, const Region &G, double eps,
                                      double coef)
{
    System_Result res;
    Vector x = x0;
    res.history.push_back(x);

    for (int k = 1; k <= MAX_ITERATIONS; k++)
    {
        Vector x_new = Phi_System(x, L, a);
        double err = coef * Diff_Norm(x_new, x);
        x = x_new;

        res.iterations = k;
        res.history.push_back(x);
        res.errors.push_back(err);

        if (!std::isfinite(err)) // расходится
            break;
        // по теореме 2.4 все приближения должны оставаться в G, иначе сходимость не гарантирована
        if (!In_Region(x, G))
        {
            res.left_region = true;
            break;
        }
        if (err <= eps)
        {
            res.converged = true;
            break;
        }
    }
    res.x = x;
    return res;
}

// Метод Ньютона для двух уравнений, формулы (2.20):
//   x1(k+1) = x1(k) - det A1 / det J
//   x2(k+1) = x2(k) - det A2 / det J
// A1 - это J, в которой первый столбец заменили на (f1, f2),
// A2 - это J, в которой второй столбец заменили на (f1, f2).
// По сути это решение системы J * dx = -f по формулам Крамера.
// Останавливаемся, когда ||x(k+1) - x(k)|| <= eps
System_Result Newton_System(const Vector &x0, double a, double eps)
{
    System_Result res;
    Vector x = x0;
    res.history.push_back(x);

    for (int k = 1; k <= MAX_ITERATIONS; k++)
    {
        Vector f = F_System(x, a);
        Matrix J = Jacobian(x, a);
        double det_J = Det2(J);
        if (std::abs(det_J) < 1e-12) // матрица Якоби вырождена, дальше считать нельзя
            break;

        Matrix A1 = {{f[0], J[0][1]}, {f[1], J[1][1]}};
        Matrix A2 = {{J[0][0], f[0]}, {J[1][0], f[1]}};
        Vector x_new = {x[0] - Det2(A1) / det_J, x[1] - Det2(A2) / det_J};

        double err = Diff_Norm(x_new, x);
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

static void Print_Matrix(const Matrix &M)
{
    for (const auto &row : M)
    {
        for (double value : row)
            std::cout << std::setw(13) << Format_Number(value);
        std::cout << "\n";
    }
}

// таблица итераций простой итерации: x(k), f(x(k)) и оценка погрешности
static void Print_Simple(const System_Result &res, double a)
{
    std::cout << std::setw(4) << "k" << std::setw(12) << "x1(k)" << std::setw(12) << "x2(k)" << std::setw(12)
              << "f1(x(k))" << std::setw(12) << "f2(x(k))" << std::setw(12) << "eps(k)" << "\n";

    for (size_t k = 0; k < res.history.size(); k++)
    {
        const Vector &x = res.history[k];
        Vector f = F_System(x, a);
        std::cout << std::setw(4) << k << std::setw(12) << Format_Number(x[0]) << std::setw(12)
                  << Format_Number(x[1]) << std::setw(12) << Format_Number(f[0], 3) << std::setw(12)
                  << Format_Number(f[1], 3);
        if (k == 0)
            std::cout << std::setw(12) << "-";
        else
            std::cout << std::setw(12) << Format_Number(res.errors[k - 1], 3);
        std::cout << "\n";
    }
    std::cout << "\n";
}

// таблица итераций Ньютона как таблица 2.5 в методичке:
// det J, det A1, det A2 считаются в точке x(k) и по ним получается x(k+1)
static void Print_Newton(const System_Result &res, double a)
{
    std::cout << std::setw(4) << "k" << std::setw(12) << "x1(k)" << std::setw(12) << "x2(k)" << std::setw(12)
              << "f1(x(k))" << std::setw(12) << "f2(x(k))" << std::setw(12) << "det J" << std::setw(12) << "det A1"
              << std::setw(12) << "det A2" << std::setw(12) << "eps(k)" << "\n";

    for (size_t k = 0; k < res.history.size(); k++)
    {
        const Vector &x = res.history[k];
        Vector f = F_System(x, a);
        Matrix J = Jacobian(x, a);
        Matrix A1 = {{f[0], J[0][1]}, {f[1], J[1][1]}};
        Matrix A2 = {{J[0][0], f[0]}, {J[1][0], f[1]}};

        std::cout << std::setw(4) << k << std::setw(12) << Format_Number(x[0]) << std::setw(12)
                  << Format_Number(x[1]) << std::setw(12) << Format_Number(f[0], 3) << std::setw(12)
                  << Format_Number(f[1], 3) << std::setw(12) << Format_Number(Det2(J)) << std::setw(12)
                  << Format_Number(Det2(A1), 3) << std::setw(12) << Format_Number(Det2(A2), 3);
        if (k == 0)
            std::cout << std::setw(12) << "-";
        else
            std::cout << std::setw(12) << Format_Number(res.errors[k - 1], 3);
        std::cout << "\n";
    }
    std::cout << "\n";
}

static void Print_Answer(const std::string &method, const System_Result &res, double a)
{
    std::cout << "Решение " << method << ": x1 = " << Format_Number(res.x[0]) << ", x2 = " << Format_Number(res.x[1])
              << ", итераций: " << res.iterations << (res.converged ? "" : " (не сошелся)") << "\n";

    // подставляем найденное решение обратно в уравнения
    Vector f = F_System(res.x, a);
    std::cout << "Проверка: f1 = " << Format_Number(f[0], 3) << ", f2 = " << Format_Number(f[1], 3) << "\n";
}

void Solve_Task_2_2(const std::string &filename)
{
    std::cout << "Задание 2.2. Метод простой итерации и метод Ньютона для системы\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле: параметр a, точность eps, начальное приближение x0 (по графику),
    // область G - границы для x1, потом границы для x2
    double a, eps;
    Vector x0(2);
    Region G;
    if (!(file >> a >> eps >> x0[0] >> x0[1] >> G.x1_min >> G.x1_max >> G.x2_min >> G.x2_max) || eps <= 0)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    std::cout << "Система (a = " << a << "):\n";
    std::cout << "f1 = a*x1^2 - x2 + x2^2 - a = 0\n";
    std::cout << "f2 = x1 - sqrt(x2 + a) + 1 = 0\n";
    std::cout << "Точность eps = " << eps << "\n";
    std::cout << "Начальное приближение (по графику): x1 = " << x0[0] << ", x2 = " << x0[1] << "\n";
    std::cout << "Область G: " << G.x1_min << " <= x1 <= " << G.x1_max << ", " << G.x2_min << " <= x2 <= "
              << G.x2_max << "\n\n";

    if (!In_Region(x0, G))
    {
        std::cout << "Начальное приближение не лежит в области G\n";
        return;
    }
    if (G.x2_min + a < 0)
    {
        std::cout << "В области G под корнем получается отрицательное число, нужно взять другую область\n";
        return;
    }

    Matrix J0 = Jacobian(x0, a);
    std::cout << "Матрица Якоби в точке x0:\n";
    Print_Matrix(J0);
    std::cout << "det J(x0) = " << Format_Number(Det2(J0)) << "\n\n";
    if (std::abs(Det2(J0)) < 1e-12)
    {
        std::cout << "Матрица Якоби вырождена, нужно взять другое начальное приближение\n";
        return;
    }

    // Метод простой итерации
    std::cout << "Метод простой итерации\n";
    std::cout << "phi(x) = x - L * f(x), где L = J^(-1)(x0):\n";
    Matrix L = Inverse2(J0);
    Print_Matrix(L);

    double q = Find_q_System(L, a, G);
    std::cout << "q = max ||phi'(x)|| по области G = " << Format_Number(q) << "\n";

    // если q < 1, погрешность оцениваем как q / (1 - q) * ||x(k) - x(k-1)||,
    // иначе просто ||x(k) - x(k-1)||
    double coef = 1.0;
    if (q < 1)
    {
        coef = q / (1 - q);
        std::cout << "q < 1, условие (2.25) выполнено\n\n";
    }
    else
    {
        std::cout << "q >= 1, сходимость не гарантирована\n\n";
    }

    System_Result simple = Simple_Iteration_System(x0, L, a, G, eps, coef);
    std::cout << "Итерации, eps(k) = " << (q < 1 ? "q / (1 - q) * " : "") << "||x(k) - x(k-1)||:\n";
    Print_Simple(simple, a);
    if (simple.left_region)
        std::cout << "Приближение вышло из области G, сходимость не гарантирована\n\n";
    else
        std::cout << "Все приближения остались в области G, значит по теореме 2.4 метод сходится\n\n";

    // Метод Ньютона
    std::cout << "Метод Ньютона\n";
    System_Result newton = Newton_System(x0, a, eps);
    std::cout << "Итерации, eps(k) = ||x(k) - x(k-1)||:\n";
    Print_Newton(newton, a);

    Print_Answer("методом простой итерации", simple, a);
    Print_Answer("методом Ньютона", newton, a);

    // как число итераций зависит от точности
    std::cout << "\nКоличество итераций при разной точности:\n";
    std::cout << std::setw(10) << "eps" << "   простая итерация   Ньютон\n";
    for (double e : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6})
    {
        System_Result s = Simple_Iteration_System(x0, L, a, G, e, coef);
        System_Result n = Newton_System(x0, a, e);
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
