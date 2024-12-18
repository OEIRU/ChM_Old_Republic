#pragma once
#include <vector>
#include "common.h"
using namespace std;
class Vector
{
public:
    Vector(int size);
    int Size();
    real EuqlideanNorm();
    real &operator()(const int index);
    friend Vector operator*(real constant, const Vector &vector);
    friend Vector operator*(const Vector &vector, real constant);
    friend real operator*(const Vector &first, const Vector &second);
    friend Vector operator+(const Vector &first, const Vector &second);

private:
    vector<real> data;
};
class Matrix
{
public:
    int Rows();
    int Columns();
    Matrix Transpose();
    Matrix(int size);
    Matrix(int rows, int columns);
    real &operator()(const int row, const int column);
    friend Matrix operator*(real constant, const Matrix &matrix);
    friend Matrix operator*(const Matrix &matrix, real constant);
    friend Vector operator*(const Matrix &matrix, Vector vector);
    friend Matrix operator*(const Matrix &first, const Matrix &second);

private:
    vector<vector<real>> data;
};