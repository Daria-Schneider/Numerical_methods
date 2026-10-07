#ifndef FILE_IO_H
#define FILE_IO_H

#include <istream>
#include "matrix.h"

bool Read_Size(std::istream &in, size_t &n);

bool Read_Matrix(std::istream &in, size_t n, Matrix &A);

bool Read_Vector(std::istream &in, size_t n, Vector &v);

#endif // FILE_IO_H
