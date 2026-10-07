#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include <string>
#include "matrix.h"

bool LU_Decomposition(const Matrix &A, Matrix &P, Matrix &L, Matrix &U, size_t &count_perm);

Vector Forward_Substitution(const Matrix &L, const Vector &b);

Vector Backward_Substitution(const Matrix &U, const Vector &z);

Vector LU_Solve(const Matrix &P, const Matrix &L, const Matrix &U, const Vector &b);

double LU_Determinant(const Matrix &U, size_t count_perm);

Matrix LU_Inverse(const Matrix &P, const Matrix &L, const Matrix &U);

void Solve_Task_1_1(const std::string &filename);

#endif // LU_DECOMPOSITION_H
