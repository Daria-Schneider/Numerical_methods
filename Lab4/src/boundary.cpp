// Задание 4.2. Краевая задача для ОДУ 2-го порядка: метод стрельбы и конечно-разностный метод,
// оценка погрешности методом Рунге-Ромберга и сравнением с точным решением
#include "boundary.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

const int MAX_ITERATIONS = 100;

// Вариант 27: (2x + 1) y'' + (2x - 1) y' - 2y = x^2 + x на отрезке [0, 1], y'(0) = 1, y'(1) + y(1) = 5.
// Делим уравнение на 2x + 1 (на [0, 1] оно больше нуля) и получаем вид (4.35):
//   y'' + p(x) y' + q(x) y = f(x)
static double p(double x)
{
    return (2 * x - 1) / (2 * x + 1);
}

static double q(double x)
{
    return -2 / (2 * x + 1);
}

static double f(double x)
{
    return (x * x + x) / (2 * x + 1);
}

// точное решение, только для проверки
static double Exact(double x)
{
    return 2 * x - 1 + std::exp(-x) + (x * x + 1) / 2;
}

// Краевые условия записаны в общем виде (4.31), условия 1-го и 2-го рода - частные случаи:
//   alpha * y'(a) + beta * y(a) = Y_A,   gamma * y'(b) + delta * y(b) = Y_B
// У нас слева y'(0) = 1 (alpha = 1, beta = 0), справа y'(1) + y(1) = 5 (gamma = delta = 1)
const double X_A = 0, X_B = 1;
const double ALPHA = 1, BETA = 0, Y_A = 1;
const double GAMMA = 1, DELTA = 1, Y_B = 5;

// Задача Коши для параметра eta, решаем методом Рунге-Кутты 4-го порядка из задания 4.1.
// Левое условие alpha y'(a) + beta y(a) = Y_A связывает y(a) и y'(a). Одно из них не известно,
// его и обозначаем eta, а второе выражаем из условия:
//  - alpha != 0 (наш случай): y(a) = eta, y'(a) = (Y_A - beta * eta) / alpha;
//  - alpha = 0 (условие 1-го рода): y(a) = Y_A / beta, y'(a) = eta, как в (4.32)
static Solution Cauchy(double eta, double h, int N)
{
    double y0, z0;
    if (ALPHA != 0)
    {
        y0 = eta;
        z0 = (Y_A - BETA * eta) / ALPHA;
    }
    else
    {
        y0 = Y_A / BETA;
        z0 = eta;
    }

    // система для z = y': y' = z, z' = f(x) - p(x) z - q(x) y
    auto rhs_y = [](double, double, double z) { return z; };
    auto rhs_z = [](double x, double y, double z) { return f(x) - p(x) * z - q(x) * y; };
    return Runge_Kutta(rhs_y, rhs_z, X_A, y0, z0, h, N);
}

// Phi(eta) из (4.33): насколько при данном eta не выполнено правое краевое условие
static double Phi(const Solution &s)
{
    return GAMMA * s.z.back() + DELTA * s.y.back() - Y_B;
}

// Метод стрельбы. Решение задачи Коши на правом конце зависит от eta, нужно подобрать eta так,
// чтобы выполнилось правое условие, т.е. решить уравнение Phi(eta) = 0.
// Производную Phi посчитать нельзя (Phi задается не формулой, а решением задачи Коши), поэтому вместо
// метода Ньютона берем метод секущих (4.34), в нем производная заменена разностью по двум последним eta:
//   eta_(j+1) = eta_j - (eta_j - eta_(j-1)) / (Phi(eta_j) - Phi(eta_(j-1))) * Phi(eta_j)
// Наше уравнение линейное, поэтому Phi(eta) - линейная функция и секущая попадает в корень за одну итерацию
bool Shooting(double h, int N, double eta0, double eta1, double eps, Solution &sol, std::vector<Shot> &shots)
{
    shots.clear();
    Solution first = Cauchy(eta0, h, N);
    shots.push_back({eta0, first.y.back(), first.z.back(), Phi(first)});
    sol = Cauchy(eta1, h, N);
    shots.push_back({eta1, sol.y.back(), sol.z.back(), Phi(sol)});

    for (int iter = 0; std::abs(shots.back().phi) >= eps; iter++)
    {
        Shot prev = shots[shots.size() - 2];
        Shot last = shots.back();
        if (iter == MAX_ITERATIONS || std::abs(last.phi - prev.phi) < 1e-15)
            return false;

        double eta = last.eta - (last.eta - prev.eta) / (last.phi - prev.phi) * last.phi;
        sol = Cauchy(eta, h, N);
        shots.push_back({eta, sol.y.back(), sol.z.back(), Phi(sol)});
    }
    return true;
}

// Метод прогонки (как в ЛР1, задание 1.2) для трехдиагональной системы
// low_k * y_(k-1) + mid_k * y_k + up_k * y_(k+1) = rhs_k
bool Sweep(const Vector &low, const Vector &mid, const Vector &up, const Vector &rhs, Vector &res)
{
    size_t n = mid.size();
    Vector P(n), Q(n);

    // прямой ход: выражаем y_k = P_k * y_(k+1) + Q_k
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

// Конечно-разностный метод. Сетка x_k = a + k*h, k = 0..N. Неизвестные y_0, ..., y_N: значения на концах
// тоже не известны, потому что в краевых условиях есть производные. Строка k системы:
//   low_k * y_(k-1) + mid_k * y_k + up_k * y_(k+1) = rhs_k
bool Finite_Difference(double h, int N, Vector &x, Vector &y, Vector &low, Vector &mid, Vector &up, Vector &rhs)
{
    x.resize(N + 1);
    for (int k = 0; k <= N; k++)
        x[k] = X_A + k * h;
    low.assign(N + 1, 0.0);
    mid.assign(N + 1, 0.0);
    up.assign(N + 1, 0.0);
    rhs.assign(N + 1, 0.0);

    // Внутренние узлы k = 1..N-1. Производные заменяем центральными разностями (4.37), это 2-й порядок:
    //   y'_k ~ (y_(k+1) - y_(k-1)) / 2h,   y''_k ~ (y_(k+1) - 2y_k + y_(k-1)) / h^2
    // Подставляем в y'' + p y' + q y = f, умножаем на h^2 и собираем коэффициенты при y_(k-1), y_k, y_(k+1),
    // получается то же, что в (4.39)
    for (int k = 1; k < N; k++)
    {
        low[k] = 1 - p(x[k]) * h / 2;
        mid[k] = -2 + h * h * q(x[k]);
        up[k] = 1 + p(x[k]) * h / 2;
        rhs[k] = h * h * f(x[k]);
    }

    // Левый конец. Если взять y'_0 ~ (y_1 - y_0) / h (4.40), будет только 1-й порядок, и из-за одной
    // граничной точки все решение станет 1-го порядка. Поэтому берем разность 2-го порядка (4.41):
    //   y'_0 ~ (-3y_0 + 4y_1 - y_2) / 2h
    // Подставляем в alpha y'_0 + beta y_0 = Y_A и умножаем на 2h:
    //   (-3 alpha + 2h beta) y_0 + 4 alpha y_1 - alpha y_2 = 2h Y_A
    // Тут три неизвестных, и матрица уже не трехдиагональная. Выражаем y_2 из уравнения для k = 1:
    //   y_2 = (rhs_1 - low_1 y_0 - mid_1 y_1) / up_1
    // подставляем, и остаются только y_0 и y_1
    mid[0] = -3 * ALPHA + 2 * h * BETA + ALPHA * low[1] / up[1];
    up[0] = 4 * ALPHA + ALPHA * mid[1] / up[1];
    rhs[0] = 2 * h * Y_A + ALPHA * rhs[1] / up[1];

    // Правый конец так же: y'_N ~ (3y_N - 4y_(N-1) + y_(N-2)) / 2h, условие gamma y'_N + delta y_N = Y_B:
    //   gamma y_(N-2) - 4 gamma y_(N-1) + (3 gamma + 2h delta) y_N = 2h Y_B
    // y_(N-2) выражаем из уравнения для k = N - 1:
    //   y_(N-2) = (rhs_(N-1) - mid_(N-1) y_(N-1) - up_(N-1) y_N) / low_(N-1)
    low[N] = -4 * GAMMA - GAMMA * mid[N - 1] / low[N - 1];
    mid[N] = 3 * GAMMA + 2 * h * DELTA - GAMMA * up[N - 1] / low[N - 1];
    rhs[N] = 2 * h * Y_B - GAMMA * rhs[N - 1] / low[N - 1];

    return Sweep(low, mid, up, rhs, y);
}

void Solve_Task_4_2(const std::string &filename)
{
    std::cout << "Задание 4.2. Краевая задача для ОДУ 2-го порядка\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    // в файле: шаг h, точность eps для метода стрельбы, два начальных значения eta0, eta1
    double h, eps, eta0, eta1;
    if (!(file >> h >> eps >> eta0 >> eta1) || h <= 0 || eps <= 0 || eta0 == eta1)
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    int N;
    if (!Steps_Count(X_A, X_B, h, N) || N < 2)
    {
        std::cout << "Шаг не укладывается на отрезке целое число раз\n";
        return;
    }
    // для Рунге-Ромберга считаем еще раз с шагом 2h; разностному методу нужно хотя бы 2 шага
    bool rr_ok = N % 2 == 0 && N >= 4;

    std::cout << "(2x + 1)y'' + (2x - 1)y' - 2y = x^2 + x, x in [0, 1]\n";
    std::cout << "y'(0) = 1, y'(1) + y(1) = 5\n";
    std::cout << "Точное решение y = 2x - 1 + e^(-x) + (x^2 + 1) / 2\n";
    std::cout << "h = " << h << ", N = " << N << "\n\n";

    // ---------- метод стрельбы ----------
    std::cout << "Метод стрельбы\n";
    std::cout << "Задача Коши y(0) = eta, y'(0) = 1 решается методом Рунге-Кутты 4-го порядка\n";
    std::cout << "Phi(eta) = y'(1) + y(1) - 5, корень ищем методом секущих, eps = " << eps << "\n";

    Solution shot;
    std::vector<Shot> shots;
    bool shot_ok = Shooting(h, N, eta0, eta1, eps, shot, shots);
    std::cout << "   j           eta          y(1)         y'(1)      Phi(eta)\n";
    for (size_t j = 0; j < shots.size(); j++)
        std::cout << std::setw(4) << j << std::setw(14) << Format_Number(shots[j].eta, 8) << std::setw(14)
                  << Format_Number(shots[j].yb, 8) << std::setw(14) << Format_Number(shots[j].zb, 8)
                  << std::setw(14) << Format_Number(shots[j].phi, 3) << "\n";

    double shot_err = 0, shot_R = 0;
    bool shot_rr = false;
    if (!shot_ok)
    {
        std::cout << "Метод секущих не сошелся\n\n";
    }
    else
    {
        std::cout << "Итераций метода секущих: " << shots.size() - 2 << "\n";
        std::cout << "Решение при eta = " << Format_Number(shots.back().eta, 8) << ":\n";
        shot_err = Print_Solution(shot.x, shot.y, Exact);

        Solution shot2;
        std::vector<Shot> shots2;
        if (rr_ok && Shooting(2 * h, N / 2, eta0, eta1, eps, shot2, shots2))
        {
            shot_rr = true;
            std::cout << "Рунге-Ромберг: решение с шагом 2h = " << 2 * h << ", p = 4, R = (y_h - y_2h) / (2^p - 1)\n";
            shot_R = Print_Runge_Romberg(shot.x, shot.y, shot2.y, 4, Exact);
        }
        std::cout << "\n";
    }

    // ---------- конечно-разностный метод ----------
    std::cout << "Конечно-разностный метод\n";
    std::cout << "y'' + p(x)y' + q(x)y = f(x), p = (2x - 1) / (2x + 1), q = -2 / (2x + 1), f = (x^2 + x) / (2x + 1)\n";
    std::cout << "Трехдиагональная система, строка k: a_k y_(k-1) + b_k y_k + c_k y_(k+1) = d_k\n";
    std::cout << "(строки 0 и N - краевые условия с разностями 2-го порядка, y_2 и y_(N-2) из них исключены)\n";

    Vector x, y, low, mid, up, rhs;
    if (!Finite_Difference(h, N, x, y, low, mid, up, rhs))
    {
        std::cout << "Прогонка не получилась (деление на ноль)\n";
        return;
    }

    // условие преобладания диагональных элементов, при нем прогонка устойчива
    bool dominant = true;
    std::cout << "   k           a_k           b_k           c_k           d_k\n";
    for (int k = 0; k <= N; k++)
    {
        std::cout << std::setw(4) << k << std::setw(14) << Format_Number(low[k]) << std::setw(14)
                  << Format_Number(mid[k]) << std::setw(14) << Format_Number(up[k]) << std::setw(14)
                  << Format_Number(rhs[k]) << "\n";
        if (std::abs(mid[k]) < std::abs(low[k]) + std::abs(up[k]))
            dominant = false;
    }
    std::cout << "Преобладание диагональных элементов |b_k| >= |a_k| + |c_k| - "
              << (dominant ? "выполнено" : "не выполнено") << "\n";

    std::cout << "Решение методом прогонки:\n";
    double fd_err = Print_Solution(x, y, Exact);

    double fd_R = 0;
    Vector x2, y2, low2, mid2, up2, rhs2;
    bool fd_rr = rr_ok && Finite_Difference(2 * h, N / 2, x2, y2, low2, mid2, up2, rhs2);
    if (fd_rr)
    {
        std::cout << "Рунге-Ромберг: решение с шагом 2h = " << 2 * h << ", p = 2, R = (y_h - y_2h) / (2^p - 1)\n";
        fd_R = Print_Runge_Romberg(x, y, y2, 2, Exact);
    }
    std::cout << "\n";

    if (!rr_ok)
        std::cout << "Число шагов нечетное или меньше 4, Рунге-Ромберг не считаем\n\n";

    std::cout << "Максимальная погрешность при h = " << h << ":\n";
    if (shot_ok)
    {
        std::cout << "  метод стрельбы: по точному решению " << Format_Number(shot_err, 3);
        if (shot_rr)
            std::cout << ", по Рунге-Ромбергу " << Format_Number(shot_R, 3);
        std::cout << "\n";
    }
    std::cout << "  конечно-разностный метод: по точному решению " << Format_Number(fd_err, 3);
    if (fd_rr)
        std::cout << ", по Рунге-Ромбергу " << Format_Number(fd_R, 3);
    std::cout << "\n";
}
