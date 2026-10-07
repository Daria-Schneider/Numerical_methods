#include "iterative_methods.h"
#include "file_io.h"
#include "lu_decomposition.h"

#include <fstream>

const int MAX_ITERATIONS = 10000;

static void Print_Iterations(const Vector &x0, const Iteration_Result &res)
{
    std::cout << std::setw(4) << "k";
    for (size_t i = 0; i < x0.size(); i++)
        std::cout << std::setw(12) << ("x" + std::to_string(i + 1));
    std::cout << std::setw(12) << "eps(k)" << "\n";

    std::cout << std::setw(4) << 0;
    for (double value : x0)
        std::cout << std::setw(12) << Format_Number(value);
    std::cout << "\n";

    for (int k = 0; k < res.iterations; k++)
    {
        std::cout << std::setw(4) << k + 1;
        for (double value : res.history[k])
            std::cout << std::setw(12) << Format_Number(value);
        std::cout << std::setw(12) << Format_Number(res.errors[k], 3) << "\n";
    }
    std::cout << "\n";
}

bool Has_Diagonal_Dominance(const Matrix &A)
{
    for (size_t i = 0; i < A.size(); i++)
    {
        double sum = 0.0;
        for (size_t j = 0; j < A.size(); j++)
            if (j != i)
                sum += std::abs(A[i][j]);
        if (std::abs(A[i][i]) <= sum)
            return false;
    }
    return true;
}

// x = beta + alpha * x, где beta_i = b_i / a_ii, alpha_ij = -a_ij / a_ii
void Build_Equivalent_System(const Matrix &A, const Vector &b, Matrix &alpha, Vector &beta)
{
    size_t n = A.size();
    alpha.assign(n, Vector(n, 0.0));
    beta.assign(n, 0.0);

    for (size_t i = 0; i < n; i++)
    {
        beta[i] = b[i] / A[i][i];
        for (size_t j = 0; j < n; j++)
            if (i != j)
                alpha[i][j] = -A[i][j] / A[i][i];
    }
}

// k + 1 >= (lg eps - lg||beta|| + lg(1 - ||alpha||)) / lg||alpha||
int A_Priori_Iterations(double norm_alpha, double norm_beta, double eps)
{
    if (norm_alpha <= 0.0 || norm_alpha >= 1.0 || norm_beta <= 0.0)
        return -1;

    double k = (std::log10(eps) - std::log10(norm_beta) + std::log10(1.0 - norm_alpha)) / std::log10(norm_alpha) - 1;
    return k < 0 ? 0 : static_cast<int>(std::ceil(k));
}

Iteration_Result Simple_Iterations(const Matrix &alpha, const Vector &beta, double eps, double coef)
{
    Iteration_Result res;
    Vector x = beta;

    for (int k = 1; k <= MAX_ITERATIONS; k++)
    {
        Vector x_new = Multiply_Vector(alpha, x);
        for (size_t i = 0; i < x.size(); i++)
            x_new[i] += beta[i];

        double err = coef * Vector_Norm(Subtract_Vector(x_new, x));
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

Iteration_Result Seidel_Method(const Matrix &alpha, const Vector &beta, double eps, double coef)
{
    Iteration_Result res;
    Vector x = beta;

    for (int k = 1; k <= MAX_ITERATIONS; k++)
    {
        // x_new меняется по ходу цикла, поэтому для j < i берутся уже новые значения
        Vector x_new = x;
        for (size_t i = 0; i < x.size(); i++)
        {
            double sum = beta[i];
            for (size_t j = 0; j < x.size(); j++)
                sum += alpha[i][j] * x_new[j];
            x_new[i] = sum;
        }

        double err = coef * Vector_Norm(Subtract_Vector(x_new, x));
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

void Solve_Task_1_3(const std::string &filename)
{
    std::cout << "Задание 1.3. Метод простых итераций и метод Зейделя\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    size_t n;
    double eps;
    Matrix A;
    Vector b;
    if (!Read_Size(file, n) || !(file >> eps) || eps <= 0 || !Read_Matrix(file, n, A) || !Read_Vector(file, n, b))
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    std::cout << "Точность eps = " << eps << "\n\n";
    std::cout << "Матрица A:\n";
    Print_Matrix(A);
    std::cout << "Вектор b:\n";
    Print_Vector(b);

    for (size_t i = 0; i < n; i++)
    {
        if (A[i][i] == 0)
        {
            std::cout << "На главной диагонали ноль, нужно переставить уравнения\n";
            return;
        }
    }

    std::cout << "\nДиагональное преобладание: " << (Has_Diagonal_Dominance(A) ? "есть" : "нет") << "\n\n";

    Matrix alpha;
    Vector beta;
    Build_Equivalent_System(A, b, alpha, beta);

    std::cout << "Приводим систему к виду x = beta + alpha * x\n";
    std::cout << "alpha:\n";
    Print_Matrix(alpha);
    std::cout << "beta:\n";
    Print_Vector(beta);

    double norm_alpha = Matrix_Norm(alpha);
    double norm_beta = Vector_Norm(beta);

    // C - верхний треугольник alpha, нужен для оценки погрешности в методе Зейделя
    Matrix C(n, Vector(n, 0.0));
    for (size_t i = 0; i < n; i++)
        for (size_t j = i; j < n; j++)
            C[i][j] = alpha[i][j];
    double norm_C = Matrix_Norm(C);

    std::cout << "\n||alpha|| = " << Format_Number(norm_alpha) << "\n";

    // если ||alpha|| < 1, то eps(k) = ||alpha|| / (1 - ||alpha||) * ||x(k) - x(k-1)||,
    // иначе просто ||x(k) - x(k-1)||
    double coef_simple = 1.0;
    double coef_seidel = 1.0;
    if (norm_alpha < 1)
    {
        coef_simple = norm_alpha / (1 - norm_alpha);
        coef_seidel = norm_C / (1 - norm_alpha);
        std::cout << "||alpha|| < 1, значит оба метода сходятся\n";
        std::cout << "Априорная оценка числа итераций: " << A_Priori_Iterations(norm_alpha, norm_beta, eps)
                  << "\n\n";
    }
    else
    {
        std::cout << "||alpha|| >= 1, сходимость не гарантирована\n\n";
    }

    Iteration_Result simple = Simple_Iterations(alpha, beta, eps, coef_simple);
    std::cout << "Метод простых итераций:\n";
    Print_Iterations(beta, simple);

    Iteration_Result seidel = Seidel_Method(alpha, beta, eps, coef_seidel);
    std::cout << "Метод Зейделя:\n";
    Print_Iterations(beta, seidel);

    std::cout << "Решение методом простых итераций (" << simple.iterations << " итераций"
              << (simple.converged ? "" : ", не сошелся") << "):\n";
    Print_Vector(simple.x);
    std::cout << "Решение методом Зейделя (" << seidel.iterations << " итераций"
              << (seidel.converged ? "" : ", не сошелся") << "):\n";
    Print_Vector(seidel.x);

    // точное решение для проверки
    Matrix P, L, U;
    size_t count_perm;
    if (LU_Decomposition(A, P, L, U, count_perm))
    {
        std::cout << "Точное решение (LU-разложение):\n";
        Print_Vector(LU_Solve(P, L, U, b));
    }

    std::cout << "\nКоличество итераций при разной точности:\n";
    std::cout << std::setw(10) << "eps" << "   простые итерации   Зейдель\n";
    for (double e : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6})
    {
        Iteration_Result s = Simple_Iterations(alpha, beta, e, coef_simple);
        Iteration_Result z = Seidel_Method(alpha, beta, e, coef_seidel);
        std::cout << std::setw(10) << e << std::setw(19) << s.iterations << std::setw(10)
                  << z.iterations << "\n";
    }
    std::cout << "\n";

    if (simple.converged && seidel.converged)
    {
        if (seidel.iterations < simple.iterations)
            std::cout << "Метод Зейделя сошелся быстрее на " << simple.iterations - seidel.iterations
                      << " итераций\n";
        else if (seidel.iterations > simple.iterations)
            std::cout << "Метод простых итераций сошелся быстрее на " << seidel.iterations - simple.iterations
                      << " итераций\n";
        else
            std::cout << "Оба метода сошлись за одинаковое число итераций\n";
    }
}
