#include <cstdlib>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

#include "boundary.h"
#include "cauchy.h"

const std::string PATH_PREFIX = LAB4_DATA_DIR;

void Run_Task(int choice)
{
    switch (choice)
    {
    case 1:
        Solve_Task_4_1(PATH_PREFIX + "input_4_1.txt");
        break;
    case 2:
        Solve_Task_4_2(PATH_PREFIX + "input_4_2.txt");
        break;
    case 0:
        std::cout << "Завершение работы программы.\n";
        break;
    default:
        std::cout << "Неверный ввод. Выберите число от 0 до 2.\n";
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

    // можно сразу передать номер задания: lab4.exe 1
    if (argc > 1)
    {
        Run_Task(std::atoi(argv[1]));
        return 0;
    }

    int choice = -1;
    while (choice != 0)
    {
        std::cout << "\nЛабораторная работа 4\n";
        std::cout << "1 - задача Коши (методы Эйлера, Рунге-Кутты и Адамса)\n";
        std::cout << "2 - краевая задача (метод стрельбы и конечно-разностный метод)\n";
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
