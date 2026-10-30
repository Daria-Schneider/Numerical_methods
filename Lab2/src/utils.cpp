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
