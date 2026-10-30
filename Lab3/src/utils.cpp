#include "utils.h"

#include <cmath>
#include <iomanip>
#include <sstream>

std::string Format_Number(double value, int precision)
{
    // погрешности округления вида 1e-17 печатаем как 0
    if (std::abs(value) < 1e-10)
        value = 0.0;

    std::ostringstream oss;
    oss << std::setprecision(precision) << value;
    return oss.str();
}

bool Read_Tokens(std::istream &in, std::vector<std::string> &tokens)
{
    tokens.clear();
    std::string line;
    while (tokens.empty() && std::getline(in, line))
    {
        std::istringstream ss(line);
        std::string word;
        while (ss >> word)
            tokens.push_back(word);
    }
    return !tokens.empty();
}

bool Read_Numbers(std::istream &in, Vector &values)
{
    std::vector<std::string> tokens;
    if (!Read_Tokens(in, tokens))
        return false;

    values.clear();
    for (const std::string &t : tokens)
    {
        std::istringstream ss(t);
        double v;
        char extra;
        if (!(ss >> v) || ss >> extra) // в слове должно быть только число
            return false;
        values.push_back(v);
    }
    return true;
}

std::string Factor_Text(double x0)
{
    if (std::abs(x0) < 1e-12)
        return "x";
    if (x0 < 0)
        return "(x + " + Format_Number(-x0) + ")";
    return "(x - " + Format_Number(x0) + ")";
}

void Add_Term(std::string &poly, double coef, const std::string &factors)
{
    if (std::abs(coef) < 1e-12)
        return;

    std::string number = Format_Number(std::abs(coef));
    if (poly.empty())
        poly = (coef < 0 ? "-" : "") + number + factors;
    else
        poly += (coef < 0 ? " - " : " + ") + number + factors;
}

std::string Equation_Text(const Vector &coefs, const std::vector<std::string> &names, double rhs)
{
    std::string s;
    for (size_t j = 0; j < coefs.size(); j++)
        Add_Term(s, coefs[j], "*" + names[j]);
    if (s.empty())
        s = "0";
    return s + " = " + Format_Number(rhs);
}
