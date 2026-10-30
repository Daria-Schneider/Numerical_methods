#ifndef UTILS_H
#define UTILS_H

#include <istream>
#include <string>
#include <vector>

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

// число для вывода в таблицы, без лишних нулей
std::string Format_Number(double value, int precision = 6);

// читает из файла следующую непустую строку и разбивает ее на слова
bool Read_Tokens(std::istream &in, std::vector<std::string> &tokens);

// читает следующую непустую строку как набор чисел
bool Read_Numbers(std::istream &in, Vector &values);

// множитель (x - x0) для записи многочлена, при x0 = 0 просто x
std::string Factor_Text(double x0);

// дописывает к записи многочлена слагаемое coef * factors, нулевые слагаемые пропускает
void Add_Term(std::string &poly, double coef, const std::string &factors);

// уравнение вида 4*c2 + 1*c3 = 1.1499 для вывода системы
std::string Equation_Text(const Vector &coefs, const std::vector<std::string> &names, double rhs);

#endif // UTILS_H
