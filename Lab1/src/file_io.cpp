#include "file_io.h"

bool Read_Size(std::istream &in, size_t &n)
{
    long long value = 0;
    if (!(in >> value) || value <= 0)
        return false;
    n = static_cast<size_t>(value);
    return true;
}

bool Read_Matrix(std::istream &in, size_t n, Matrix &A)
{
    A.assign(n, Vector(n, 0.0));
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (!(in >> A[i][j]))
                return false;
    return true;
}

bool Read_Vector(std::istream &in, size_t n, Vector &v)
{
    v.assign(n, 0.0);
    for (size_t i = 0; i < n; i++)
        if (!(in >> v[i]))
            return false;
    return true;
}
