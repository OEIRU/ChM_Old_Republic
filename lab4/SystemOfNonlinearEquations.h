#pragma once
#include "LinearAlgebra.h"

class VectorOfFunctions
{
public:
    virtual Vector ComputeInPoint(Vector point) = 0;
};
class DisjointCircles : public VectorOfFunctions
{
private:
    int size = 2;

public:
    Vector ComputeInPoint(Vector point) override;
};
class IntersectingCirclesAtPoint : public VectorOfFunctions
{
private:
    int size = 2;

public:
    Vector ComputeInPoint(Vector point) override;
};
class IntersectingCircles : public VectorOfFunctions
{
private:
    int size = 2;

public:
    Vector ComputeInPoint(Vector point) override;
};
struct SystemParameters
{
    int n;
    int m;
    int maxiter;
    int maxiterBeta;
    double epsF;
    double epsBeta;
    VectorOfFunctions *F;
    Vector x0 = Vector(1);
    SystemParameters(int n, int m, int maxiter, int maxiterBeta, double epsF,
        double epsBeta, VectorOfFunctions *Function, Vector x0);
};
class Squaring
{
public:
    virtual void LeadToSquare(Matrix &matrix, Vector &vector);
};
class ExcludingRows : public Squaring
{
public:
    void LeadToSquare(Matrix &matrix, Vector &vector) override;
};
class Convolution : public Squaring
{
public:
    void LeadToSquare(Matrix &matrix, Vector &vector) override;
};
class SystemOfNonlinearEquations
{
private:
    int n;
    int m;
    int maxiter;
    int maxiterBeta;
    double epsF;
    double epsBeta;
    VectorOfFunctions *F;
    Vector x0 = Vector(1);
    Squaring *squaring;

public:
    SystemOfNonlinearEquations(struct SystemParameters parameters, Squaring *
                                                                       squaring);
    Vector Solve();

private:
    Matrix FormJacobiMatrix(Vector x);
    Vector ComputeDirectionByGauss(Matrix matrix, Vector vector);
};