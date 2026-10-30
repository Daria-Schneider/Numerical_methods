#include <cstdlib>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

#include "differentiation.h"
#include "integration.h"
#include "interpolation.h"
#include "least_squares.h"
#include "spline.h"

const std::string PATH_PREFIX = LAB3_DATA_DIR;

void Run_Task(int choice)
{
    switch (choice)
    {
    case 1:
        Solve_Task_3_1(PATH_PREFIX + "input_3_1.txt");
        break;
    case 2:
        Solve_Task_3_2(PATH_PREFIX + "input_3_2.txt");
        break;
    case 3:
        Solve_Task_3_3(PATH_PREFIX + "input_3_3.txt");
        break;
    case 4:
        Solve_Task_3_4(PATH_PREFIX + "input_3_4.txt");
        break;
    case 5:
        Solve_Task_3_5(PATH_PREFIX + "input_3_5.txt");
        break;
    case 0:
        std::cout << "Завершение работы программы.\n";
        break;
    default:
        std::cout << "Неверный ввод. Выберите число от 0 до 5.\n";
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

    // можно сразу передать номер задания: lab3.exe 1
    if (argc > 1)
    {
        Run_Task(std::atoi(argv[1]));
        return 0;
    }

    int choice = -1;
    while (choice != 0)
    {
        std::cout << "\nЛабораторная работа 3\n";
        std::cout << "1 - интерполяционные многочлены Лагранжа и Ньютона\n";
        std::cout << "2 - кубический сплайн\n";
        std::cout << "3 - метод наименьших квадратов\n";
        std::cout << "4 - численное дифференцирование\n";
        std::cout << "5 - численное интегрирование\n";
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
