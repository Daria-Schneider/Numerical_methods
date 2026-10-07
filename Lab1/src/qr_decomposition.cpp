#include "qr_decomposition.h"
#include "file_io.h"

#include <fstream>

using Complex = std::complex<double>;

const int MAX_ITERATIONS = 1000;

// норма элементов столбца col, начиная со строки from
static double Column_Norm(const Matrix &A, size_t col, size_t from)
{
    double sum = 0.0;
    for (size_t l = from; l < A.size(); l++)
        sum += A[l][col] * A[l][col];
    return std::sqrt(sum);
}

// собственные значения блока 2x2 в строках/столбцах j, j+1:
// (a_jj - l)(a_j+1,j+1 - l) = a_j,j+1 * a_j+1,j
// возвращает true, если корни комплексные
static bool Block_Eigenvalues(const Matrix &A, size_t j, Complex &l1, Complex &l2)
{
    double a = A[j][j], b = A[j][j + 1];
    double c = A[j + 1][j], d = A[j + 1][j + 1];

    double D = (a - d) * (a - d) + 4 * b * c;
    if (D >= 0)
    {
        l1 = (a + d + std::sqrt(D)) / 2;
        l2 = (a + d - std::sqrt(D)) / 2;
        return false;
    }
    l1 = Complex((a + d) / 2, std::sqrt(-D) / 2);
    l2 = Complex((a + d) / 2, -std::sqrt(-D) / 2);
    return true;
}

// QR-разложение через матрицы Хаусхолдера H = E - 2 * v * v^T / (v^T * v)
void QR_Decomposition(const Matrix &A, Matrix &Q, Matrix &R)
{
    size_t n = A.size();
    Q = Identity_Matrix(n);
    R = A;

    for (size_t k = 0; k + 1 < n; k++)
    {
        double norm = Column_Norm(R, k, k);
        if (norm < 1e-15)
            continue;

        // v_k = a_kk + sign(a_kk) * ||a||, ниже - элементы столбца, выше - нули
        Vector v(n, 0.0);
        v[k] = R[k][k] + (R[k][k] >= 0 ? norm : -norm);
        for (size_t i = k + 1; i < n; i++)
            v[i] = R[i][k];

        double vv = 0.0;
        for (double value : v)
            vv += value * value;

        Matrix H = Identity_Matrix(n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                H[i][j] -= 2 * v[i] * v[j] / vv;

        R = Multiply_Matrix(H, R);
        Q = Multiply_Matrix(Q, H);

        for (size_t i = k + 1; i < n; i++)
            R[i][k] = 0.0; // там остаются числа порядка 1e-16
    }
}

QR_Result QR_Algorithm(const Matrix &A, double eps, bool print)
{
    size_t n = A.size();
    QR_Result res;
    Matrix Ak = A;
    std::vector<Complex> prev;

    if (print)
    {
        // заголовок пишем без setw: он считает байты, а русская буква занимает 2 байта
        std::cout << std::setw(4) << "k";
        for (size_t j = 0; j + 1 < n; j++)
            std::cout << "   столбец " << j + 1;
        std::cout << "   собственные значения\n";
    }

    for (int iter = 1; iter <= MAX_ITERATIONS; iter++)
    {
        Matrix Q, R;
        QR_Decomposition(Ak, Q, R);
        Ak = Multiply_Matrix(R, Q);
        res.iterations = iter;

        // смотрим, где на диагонали вещественное значение, а где блок 2x2 с комплексной парой
        std::vector<Complex> lambda(n);
        std::vector<bool> is_pair(n, false);
        for (size_t j = 0; j < n;)
        {
            Complex l1, l2;
            if (j + 1 < n && Column_Norm(Ak, j, j + 1) > eps && Block_Eigenvalues(Ak, j, l1, l2))
            {
                lambda[j] = l1;
                lambda[j + 1] = l2;
                is_pair[j] = true;
                j += 2;
            }
            else
            {
                lambda[j] = Ak[j][j];
                j += 1;
            }
        }

        // проверка окончания:
        // для вещественного значения - элементы под диагональю в столбце <= eps,
        // для комплексной пары - |lambda(k) - lambda(k-1)| <= eps
        bool done = true;
        for (size_t j = 0; j < n;)
        {
            if (is_pair[j])
            {
                if (Column_Norm(Ak, j, j + 2) > eps || Column_Norm(Ak, j + 1, j + 2) > eps)
                    done = false;
                if (prev.empty() || std::abs(lambda[j] - prev[j]) > eps || std::abs(lambda[j + 1] - prev[j + 1]) > eps)
                    done = false;
                j += 2;
            }
            else
            {
                if (Column_Norm(Ak, j, j + 1) > eps)
                    done = false;
                j += 1;
            }
        }

        if (print)
        {
            std::cout << std::setw(4) << iter;
            for (size_t j = 0; j + 1 < n; j++)
                std::cout << std::setw(12) << Format_Number(Column_Norm(Ak, j, j + 1), 3);
            std::cout << "   ";
            for (size_t j = 0; j < n; j++)
                std::cout << Format_Complex(lambda[j], 4) << (j + 1 < n ? ";  " : "\n");
        }

        prev = lambda;
        if (done)
        {
            res.converged = true;
            break;
        }
    }

    res.A_final = Ak;
    res.eigenvalues = prev;
    return res;
}

void Solve_Task_1_5(const std::string &filename)
{
    std::cout << "Задание 1.5. QR-алгоритм\n\n";

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

    Matrix Q, R;
    QR_Decomposition(A, Q, R);
    std::cout << "Матрица Q:\n";
    Print_Matrix(Q);
    std::cout << "Матрица R:\n";
    Print_Matrix(R);

    Matrix QR = Multiply_Matrix(Q, R);
    std::cout << "Q * R:\n";
    Print_Matrix(QR);
    std::cout << "Q * R = A: " << (Equal_Matrix(QR, A) ? "верно" : "неверно") << "\n\n";

    std::cout << "Итерации (для каждого столбца - норма элементов под диагональю):\n";
    QR_Result res = QR_Algorithm(A, eps, true);

    std::cout << "\nКоличество итераций: " << res.iterations << (res.converged ? "" : " (точность не достигнута)")
              << "\n\n";
    std::cout << "A(" << res.iterations << "):\n";
    Print_Matrix(res.A_final);

    std::cout << "Собственные значения:\n";
    for (size_t i = 0; i < n; i++)
        std::cout << "  lambda" << i + 1 << " = " << Format_Complex(res.eigenvalues[i]) << "\n";
}
