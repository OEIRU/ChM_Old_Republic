#ifndef MAIN_H
#define MAIN_H

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#define EPS 1E-14  // Точность сравнения
#define _div(dst, a, b) { if (((b) > (T)0 ? (b) : -(b)) < (T)EPS) return false; dst = (a) / (b); }
#define elementsInLine(i) (a.ia[i + 1] - a.ia[i])


// Абсолютное значение для double
inline double dabs(double value) { return value > 0 ? value : -value; }


// Матрица в профильном формате
template<typename T>
struct SkylineStorageMatrix
{
    std::vector<int> ia;  // Индексы начала строк
    std::vector<T> al;    // Нижняя часть матрицы (левая)
    std::vector<T> au;    // Верхняя часть матрицы (правая)
    std::vector<T> di;    // Диагональные элементы
};

// Генерация матрицы Гильберта
template<typename T>
void generateHilbert(SkylineStorageMatrix<T>& dst, std::vector<T>& dstr, int n);

// Тесты для матрицы Гильберта
template<typename T>
void hilbertTestSeries(int n1, int n2);

// Ввод/вывод для матрицы в профильном формате
template<typename T>
void readSkyline(std::istream& from, SkylineStorageMatrix<T>& to);

template<typename T>
void skylineToSquare(SkylineStorageMatrix<T>& a, std::vector<std::vector<T>>& s);

// Базовая работа с векторами и матрицами
template<typename T>
void readVector(std::istream& from, std::vector<T>& to);

template<typename T>
void printMatrix(const std::vector<std::vector<T>>& a);

template<typename T, typename T2=T>
bool LusqSkyline(SkylineStorageMatrix<T>& a);

template<typename T>
void printSkyline(SkylineStorageMatrix<T>& a);

template<typename T>
void writeVector(std::vector<T>& from, std::ostream& to);

template<typename T>
void writeSkyline(SkylineStorageMatrix<T>& from, std::ostream& to);

// LU разложение и решение системы
template<typename T, typename T2=T>
bool findY(std::vector<T>& y, SkylineStorageMatrix<T>& a, std::vector<T>& b);

template<typename T, typename T2=T>
bool findX(std::vector<T>& x, SkylineStorageMatrix<T>& a, std::vector<T>& y);

template<typename T, typename T2=T>
bool solve(std::vector<T>& x, SkylineStorageMatrix<T>& a, std::vector<T>& b);

// Решение методом Гаусса
template<typename T>
bool solveGauss(std::vector<std::vector<T>> A, std::vector<T>& x, std::vector<T> B);

// Тестирование
void compareGaussLusq(int tests);
void conditionNumberTestSeries(int tests);

// Меню
void showMenu();

#endif // MAIN_H
