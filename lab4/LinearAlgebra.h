#pragma once
#include <vector>
using namespace std;
class Vector
{
public:
    Vector(int size);
    int Size();
    double EuqlideanNorm();
    double &operator()(const int index);
    friend Vector operator*(double constant, const Vector &vector);
    friend Vector operator*(const Vector &vector, double constant);
    friend double operator*(const Vector &first, const Vector &second);
    friend Vector operator+(const Vector &first, const Vector &second);

private:
    vector<double> data;
};
class Matrix
{
public:
    int Rows();
    int Columns();
    Matrix Transpose();
    Matrix(int size);
    Matrix(int rows, int columns);
    double &operator()(const int row, const int column);
    friend Matrix operator*(double constant, const Matrix &matrix);
    friend Matrix operator*(const Matrix &matrix, double constant);
    friend Vector operator*(const Matrix &matrix, Vector vector);
    friend Matrix operator*(const Matrix &first, const Matrix &second);

private:
    vector<vector<double>> data;
};