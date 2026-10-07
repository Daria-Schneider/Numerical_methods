#include "jacobi_eigen.h"
#include "file_io.h"

#include <fstream>

const int MAX_ITERATIONS = 10000;

bool Is_Symmetric(const Matrix &A)
{
    for (size_t i = 0; i < A.size(); i++)
        for (size_t j = i + 1; j < A.size(); j++)
            if (std::abs(A[i][j] - A[j][i]) > 1e-9)
                return false;
    return true;
}

// t(A) - корень из суммы квадратов элементов над главной диагональю
double Off_Diagonal_Norm(const Matrix &A)
{
    double sum = 0.0;
    for (size_t l = 0; l < A.size(); l++)
        for (size_t m = l + 1; m < A.size(); m++)
            sum += A[l][m] * A[l][m];
    return std::sqrt(sum);
}

Jacobi_Result Jacobi_Rotation(const Matrix &A, double eps)
{
    const double PI = std::acos(-1.0);
    size_t n = A.size();

    Jacobi_Result res;
    Matrix Ak = A;
    Matrix V = Identity_Matrix(n);

    double t = Off_Diagonal_Norm(Ak);
    res.t0 = t;

    while (t > eps && res.iterations < MAX_ITERATIONS)
    {
        // максимальный по модулю элемент над диагональю
        size_t i = 0, j = 1;
        for (size_t l = 0; l < n; l++)
            for (size_t m = l + 1; m < n; m++)
                if (std::abs(Ak[l][m]) > std::abs(Ak[i][j]))
                {
                    i = l;
                    j = m;
                }

        // угол поворота
        double phi;
        if (Ak[i][i] == Ak[j][j])
            phi = PI / 4;
        else
            phi = 0.5 * std::atan(2 * Ak[i][j] / (Ak[i][i] - Ak[j][j]));

        // матрица вращения
        Matrix U = Identity_Matrix(n);
        U[i][i] = std::cos(phi);
        U[i][j] = -std::sin(phi);
        U[j][i] = std::sin(phi);
        U[j][j] = std::cos(phi);

        // A(k+1) = U^T * A(k) * U
        Ak = Multiply_Matrix(Multiply_Matrix(Transpose(U), Ak), U);
        V = Multiply_Matrix(V, U);

        res.iterations++;
        t = Off_Diagonal_Norm(Ak);
        res.steps.push_back({res.iterations, i, j, phi, t});
    }

    res.converged = (t <= eps);
    res.A_final = Ak;
    res.eigenvectors = V;
    for (size_t i = 0; i < n; i++)
        res.eigenvalues.push_back(Ak[i][i]);
    return res;
}

void Solve_Task_1_4(const std::string &filename)
{
    std::cout << "Задание 1.4. Метод вращений Якоби\n\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    size_t n;
    double eps;
    Matrix A;
    if (!Read_Size(file, n) || !(file >> eps) || eps <= 0 || !Read_Matrix(file, n, A))
    {
        std::cerr << "Ошибка: неверные данные в файле " << filename << "\n";
        return;
    }

    std::cout << "Точность eps = " << eps << "\n\n";
    std::cout << "Матрица A:\n";
    Print_Matrix(A);

    if (!Is_Symmetric(A))
    {
        std::cout << "Матрица не симметрическая, метод вращений не применим\n";
        return;
    }

    Jacobi_Result res = Jacobi_Rotation(A, eps);

    std::cout << "Итерации:\n";
    std::cout << std::setw(4) << "k" << std::setw(8) << "(i,j)" << std::setw(12) << "phi" << std::setw(12) << "t(A)"
              << "\n";
    std::cout << std::setw(4) << 0 << std::setw(8) << "-" << std::setw(12) << "-" << std::setw(12)
              << Format_Number(res.t0, 4) << "\n";
    for (const auto &s : res.steps)
    {
        std::string pos = "(" + std::to_string(s.i + 1) + "," + std::to_string(s.j + 1) + ")";
        std::cout << std::setw(4) << s.k << std::setw(8) << pos << std::setw(12) << Format_Number(s.phi)
                  << std::setw(12) << Format_Number(s.t, 4) << "\n";
    }
    std::cout << "\nКоличество итераций: " << res.iterations << (res.converged ? "" : " (точность не достигнута)")
              << "\n\n";

    std::cout << "Матрица после вращений:\n";
    Print_Matrix(res.A_final);

    std::cout << "Собственные значения:\n";
    for (size_t i = 0; i < n; i++)
        std::cout << "  lambda" << i + 1 << " = " << Format_Number(res.eigenvalues[i]) << "\n";

    std::cout << "\nСобственные векторы:\n";
    for (size_t j = 0; j < n; j++)
    {
        std::cout << "  x" << j + 1 << " = (";
        for (size_t i = 0; i < n; i++)
            std::cout << Format_Number(res.eigenvectors[i][j]) << (i + 1 < n ? ", " : ")\n");
    }

    std::cout << "\nМатрица собственных векторов U:\n";
    Print_Matrix(res.eigenvectors);

    // проверка A * U = U * Lambda
    Matrix Lambda(n, Vector(n, 0.0));
    for (size_t i = 0; i < n; i++)
        Lambda[i][i] = res.eigenvalues[i];

    Matrix AU = Multiply_Matrix(A, res.eigenvectors);
    Matrix ULambda = Multiply_Matrix(res.eigenvectors, Lambda);
    std::cout << "A * U:\n";
    Print_Matrix(AU);
    std::cout << "U * Lambda:\n";
    Print_Matrix(ULambda);

    // разность A*U - U*Lambda по модулю не больше sqrt(2) * t(A), а t(A) <= eps
    bool ok = Equal_Matrix(AU, ULambda, std::sqrt(2.0) * eps);
    std::cout << "A * U = U * Lambda: " << (ok ? "верно" : "неверно") << "\n\n";

    std::cout << "Скалярные произведения собственных векторов:\n";
    for (size_t l = 0; l < n; l++)
        for (size_t m = l + 1; m < n; m++)
        {
            double dot = 0.0;
            for (size_t i = 0; i < n; i++)
                dot += res.eigenvectors[i][l] * res.eigenvectors[i][m];
            std::cout << "  (x" << l + 1 << ", x" << m + 1 << ") = " << Format_Number(dot) << "\n";
        }

    std::cout << "\nКоличество итераций при разной точности:\n";
    std::cout << std::setw(10) << "eps" << "   итераций\n";
    for (double e : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10})
        std::cout << std::setw(10) << e << std::setw(11) << Jacobi_Rotation(A, e).iterations << "\n";
}
