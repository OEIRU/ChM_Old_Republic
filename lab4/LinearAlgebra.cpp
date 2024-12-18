#include "LinearAlgebra.h"
#include "common.h"
#include <fstream>
Vector::Vector(int size)
{
    data.resize(size, 0.0);
}
int Vector::Size()
{
    return data.size();
}
double Vector::EuqlideanNorm()
{
    double sum = 0.0;
    for (int i = 0; i < data.size(); i++)
    {
        sum += data[i] * data[i];
    }
    return sqrt(sum);
}
double&Vector::operator()(int index)
{
    return data[index];
}
Vector operator*(double constant, const Vector &vector)
{
    Vector result = vector;
    for (int i = 0; i < result.Size(); i++)
    {
        result.data[i] *= constant;
    }
    return result;
}
Vector operator*(const Vector &vector, double constant)
{
    Vector result = vector;
    for (int i = 0; i < result.Size(); i++)
    {
        result(i) *= constant;
    }
    return result;
}
double operator*(const Vector &first, const Vector &second)
{
    double result = 0.0;
    for (int i = 0; i < second.data.size(); i++)
    {
        result += first.data[i] * second.data[i];
    }
    return result;
}
Vector operator+(const Vector &first, const Vector &second)
{
    Vector result(first.data.size());
    for (int i = 0; i < first.data.size(); i++)
    {
        result.data[i] = first.data[i] + second.data[i];
    }
    return result;
}
Matrix::Matrix(int size)
{
    data.resize(size);
    for (int i = 0; i < size; i++)
    {
        data[i].resize(size, 0.0);
    }
}
Matrix::Matrix(int rows, int colums)
{
    data.resize(rows);
    for (int i = 0; i < rows; i++)
    {
        data[i].resize(colums, 0.0);
    }
}
int Matrix::Columns()
{
    return data[0].size();
}
int Matrix::Rows()
{
    return data.size();
}
Matrix Matrix::Transpose()
{
    Matrix result(Columns(), Rows());
    for (int i = 0; i < Rows(); i++)
    {
        for (int j = 0; j < Columns(); j++)
        {
            result.data[j][i] = data[i][j];
        }
    }
    return result;
}
Matrix operator*(double constant, const Matrix &matrix)
{
    Matrix result(matrix.data.size());
    for (int i = 0; i < result.data.size(); i++)
    {
        for (int j = 0; j < result.data[i].size(); j++)
        {
            result.data[i][j] *= constant;
        }
    }
    return result;
}
Matrix operator*(const Matrix &matrix, double constant)
{
    Matrix result(matrix.data.size());
    for (int i = 0; i < result.data.size(); i++)
    {
        for (int j = 0; j < result.data[i].size(); j++)
        {
            result.data[i][j] *= constant;
        }
    }
    return result;
}
Matrix operator*(const Matrix &first, const Matrix &second)
{
    Matrix result(first.data.size());
    for (int i = 0; i < first.data.size(); i++)
    {
        for (int j = 0; j < first.data.size(); j++)
        {
            for (int k = 0; k < second.data.size(); k++)
            {
                result.data[i][j] += first.data[i][k] *
                                     second.data[k][j];
            }
        }
    }
    return result;
}
double &Matrix::operator()(const int row, const int column)
{
    return data[row][column];
}
Vector operator*(const Matrix &matrix, Vector vector)
{
    Vector result(matrix.data.size());
    for (int i = 0; i < matrix.data.size(); i++)
    {
        for (int j = 0; j < matrix.data[0].size(); j++)
        {
            result(i) += matrix.data[i][j] * vector(j);
        }
    }
    return result;
}