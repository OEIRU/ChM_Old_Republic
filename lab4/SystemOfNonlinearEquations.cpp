#pragma once
#include "SystemOfNonlinearEquations.h"
#include <iostream>
SystemParameters::SystemParameters(int n, int m, int maxiter, int maxiterBeta,
                                   real epsF, real epsBeta, VectorOfFunctions *Function, Vector x0)
{
    this->n = n;
    this->m = m;
    this->maxiter = maxiter;
    this->maxiterBeta = maxiterBeta;
    this->epsF = epsF;
    this->epsBeta = epsBeta;
    this->F = Function;
    this->x0 = x0;
}
Vector DisjointCircles::ComputeInPoint(Vector point)
{
    Vector value(size);
    value(0) = pow(point(0) + 2, 2) + pow(point(1) - 2, 2) - 4;
    value(1) = pow(point(0) - 2, 2) + pow(point(1) - 2, 2) - 4;
    return value;
}
Vector IntersectingCirclesAtPoint::ComputeInPoint(Vector point)
{
    return Vector(2);
}
Vector IntersectingCircles::ComputeInPoint(Vector point)
{
    return Vector(2);
}
void Squaring::LeadToSquare(Matrix &matrix, Vector &vector)
{
    Matrix temp = matrix.Transpose();
    matrix = temp * matrix;
    vector = (-1.0 * temp) * vector;
}
void ExcludingRows::LeadToSquare(Matrix &matrix, Vector &vector)
{
    int m = matrix.Rows();
    int n = matrix.Columns();
    int rowsToDelete = m - n;
    int row;
    real min;
    for (int i = 0; i < rowsToDelete; i++)
    {
        row = 0;
        min = fabs(vector(0));
        for (int j = 0; j < m - i; j++)
        {
            if (fabs(vector(j)) < min)
            {
                min = fabs(vector(j));
                row = j;
            }
        }
        swap(vector(row), vector(m - i - 1));
        for (int j = 0; j < n; j++)
        {
            swap(matrix(row, j), matrix(m - i - 1, j));
        }
    }
}
void Convolution::LeadToSquare(Matrix &matrix, Vector &vector)
{
    int m = matrix.Rows();
    int n = matrix.Columns();
    int rowsToDelete = m - n + 1;
    int row;
    real min;
    real sum = 0.0;
    Vector rowsSum(n);
    for (int i = 0; i < rowsToDelete; i++)
    {
        row = 0;
        min = fabs(vector(0));
        for (int j = 0; j < m - i; j++)
        {
            if (fabs(vector(j)) < min)
            {
                min = fabs(vector(j));
                row = j;
            }
        }
        sum += pow(vector(row), 2);
        for (int j = 0; j < n; j++)
        {
            rowsSum(j) += pow(matrix(row, j), 2);
        }
        swap(vector(row), vector(m - i - 1));
        for (int j = 0; j < n; j++)
        {
            swap(matrix(row, j), matrix(m - i - 1, j));
        }
    }
    vector(n - 1) = sum;
    for (int i = 0; i < n; i++)
    {
        matrix(n - 1, i) = rowsSum(i);
    }
}
SystemOfNonlinearEquations::SystemOfNonlinearEquations(struct SystemParameters
                                                           parameters,
                                                       Squaring *squaring)
{
    this->n = parameters.n;
    this->m = parameters.m;
    this->maxiter = parameters.maxiter;
    this->maxiterBeta = parameters.maxiterBeta;
    this->epsF = parameters.epsF;
    this->epsBeta = parameters.epsBeta;
    this->F = parameters.F;
    this->x0 = parameters.x0;
    this->squaring = squaring;
}
Matrix SystemOfNonlinearEquations::FormJacobiMatrix(Vector x)
{
    Matrix Jacobi(m, n);
    real h = 1e-10;
    Vector temp = x;
    Vector Fp = F->ComputeInPoint(x);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp(j) += h;
            Vector Fp_h = F->ComputeInPoint(temp);
            for (int k = 0; k < m; k++)
            {
                Fp_h(k) -= Fp(k);
            }
            Jacobi(i, j) = Fp_h(i) / h;
            temp(j) = x(j);
        }
    }
    return Jacobi;
}
Vector SystemOfNonlinearEquations::ComputeDirectionByGauss(Matrix matrix, Vector
                                                                              vector)
{
    vector = -1 * vector;
    int i;
    for (i = 0; i < n; i++)
    {
        real mainElement = 0.0;
        int row = 0;
        for (int j = i; j < n; j++)
        {
            if (mainElement < fabs(matrix(j, i)))
            {
                mainElement = matrix(j, i);
                row = j;
            }
        }
        if (row != i)
        {
            swap(vector(i), vector(row));
            for (int j = 0; j < n; j++)
            {
                swap(matrix(i, j), matrix(row, j));
            }
        }
        vector(i) /= mainElement;
        for (int j = i + 1; j < n; j++)
        {
            matrix(i, j) /= mainElement;
        }
        for (int j = i + 1; j < n; j++)
        {
            mainElement = matrix(j, i);
            for (int k = i; k < n; k++)
            {
                matrix(j, k) -= mainElement * matrix(i, k);
            }
            vector(j) -= mainElement * vector(i);
        }
    }
    for (i -= 2; i >= 0; i--)
    {
        for (int j = i + 1; j < n; j++)
        {
            vector(i) -= vector(j) * matrix(i, j);
        }
    }
    return vector;
}
Vector SystemOfNonlinearEquations::Solve()
{
    real F0Norm = F->ComputeInPoint(x0).EuqlideanNorm();
    real FNorm;
    real FkNorm = F0Norm;
    Vector xk = x0;
    Vector xk1(n);
    real discrepancy = FkNorm / F0Norm;
    for (int k = 0; k < maxiter && discrepancy > epsF; k++)
    {
        Vector Fk = F->ComputeInPoint(xk);
        Matrix Jacobi = FormJacobiMatrix(xk);
        if (m != n)
        {
            squaring->LeadToSquare(Jacobi, Fk);
        }
        Vector dx = ComputeDirectionByGauss(Jacobi, Fk);
        real beta = 1.0;
        FNorm = FkNorm;
        for (int v = 0; v < maxiterBeta && beta > epsBeta; v++)
        {
            for (int i = 0; i < n; i++)
            {
                xk1(i) = xk(i) + beta * dx(i);
            }
            FkNorm = F->ComputeInPoint(xk1).EuqlideanNorm();
            if (FkNorm < FNorm)
            {
                break;
            }
            beta /= 2;
        }
        for (int i = 0; i < n; i++)
        {
            xk(i) = xk1(i);
        }
        discrepancy = FkNorm / F0Norm;
        cout << "beta: " << beta << endl;
        cout << "discrepancy: " << discrepancy << endl;
        cout << "xk: ";
        for (int i = 0; i < n; i++)
        {
            cout << fixed << xk(i) << " ";
        }
        cout << endl
             << endl;
    }
    return xk;
}