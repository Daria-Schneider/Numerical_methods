#include <cstdlib>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

#include "lu_decomposition.h"
#include "tridiagonal_sweep.h"
#include "iterative_methods.h"
#include "jacobi_eigen.h"
#include "qr_decomposition.h"

const std::string PATH_PREFIX = LAB1_DATA_DIR;

void Run_Task(int choice)
{
    switch (choice)
    {
    case 1:
        Solve_Task_1_1(PATH_PREFIX + "input_1_1.txt");
        break;
    case 2:
        Solve_Task_1_2(PATH_PREFIX + "input_1_2.txt");
        break;
    case 3:
        Solve_Task_1_3(PATH_PREFIX + "input_1_3.txt");
        break;
    case 4:
        Solve_Task_1_4(PATH_PREFIX + "input_1_4.txt");
        break;
    case 5:
        Solve_Task_1_5(PATH_PREFIX + "input_1_5.txt");
        break;
    case 6:
        Solve_Task_1_5(PATH_PREFIX + "input_1_5_complex.txt");
        break;
    case 0:
        std::cout << "Завершение работы программы.\n";
        break;
    default:
        std::cout << "Неверный ввод. Выберите число от 0 до 6.\n";
        break;
    }
}

int main(int argc, char *argv[])
{
#ifdef _WIN32
    // чтобы в консоли нормально отображался русский текст
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    // можно сразу передать номер задания: lab1.exe 3
    if (argc > 1)
    {
        Run_Task(std::atoi(argv[1]));
        return 0;
    }

    int choice = -1;
    while (choice != 0)
    {
        std::cout << "\nЛабораторная работа 1\n";
        std::cout << "1 - LU-разложение\n";
        std::cout << "2 - метод прогонки\n";
        std::cout << "3 - метод простых итераций и метод Зейделя\n";
        std::cout << "4 - метод вращений\n";
        std::cout << "5 - QR-алгоритм\n";
        std::cout << "6 - QR-алгоритм, пример с комплексными корнями\n";
        std::cout << "0 - выход\n";
        std::cout << "Номер задания: ";

        if (!(std::cin >> choice))
        {
            if (std::cin.eof())
                break;
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            choice = -1;
        }

        std::cout << "\n";
        Run_Task(choice);
    }

    return 0;
}
