#include "LinearAlgebra.h"
#include "SystemOfNonlinearEquations.h"
#include <iostream>

int main() {
    // ��������� ������
    int n = 2; // ����������� �������
    int m = 2; // ���������� ���������
    int maxIter = 100; // ������������ ���������� ��������
    int maxIterBeta = 50; // ������������ ���������� �������� ��� ��������� ����
    double epsF = 1e-6; // �������� �� ����� �������
    double epsBeta = 1e-6; // �������� �� ����

    // ������ ��������� �����������, ��� � ��������
    Vector x0(n);
    x0(0) = -3; // x1 ���������
    x0(1) = 0;  // x2 ���������

    // ���������� ������� ���������� ��������� (�������������� ����������)
    DisjointCircles function; // ������������ ����� ��� �������� �����������

    // ������� ��������� �������
    SystemParameters params(n, m, maxIter, maxIterBeta, epsF, epsBeta, &function, x0);

    // �������� ����� ��������� (�������� Symmetrization �� Convolution ��� ExcludingRows)
    Convolution squaring; // ������������ ����� �������

    // ������� ������ ������� �������
    SystemOfNonlinearEquations solver(params, &squaring);

    // �������� ������ �������
    Vector solution = solver.Solve();

    // ������� ���������
    std::cout << "������� �������:" << std::endl;
    std::cout << "x1 = " << solution(0) << std::endl;
    std::cout << "x2 = " << solution(1) << std::endl;

    // ����������� ��������� (��������, ����� ����� ������� � �������� �����)
    double finalDiscrepancy = function.ComputeInPoint(solution).EuqlideanNorm();
    std::cout << "����� ������� � �������: " << finalDiscrepancy << std::endl;

    return 0;
}
