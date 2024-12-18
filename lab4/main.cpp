#include "LinearAlgebra.h"
#include "SystemOfNonlinearEquations.h"
#include <iostream>

int main() {
    // Параметры задачи
    int n = 2; // Размерность системы
    int m = 2; // Количество уравнений
    int maxIter = 100; // Максимальное количество итераций
    int maxIterBeta = 50; // Максимальное количество итераций для параметра бета
    double epsF = 1e-6; // Точность по норме функции
    double epsBeta = 1e-6; // Точность по бета

    // Задаем начальное приближение, как в графиках
    Vector x0(n);
    x0(0) = -3; // x1 начальное
    x0(1) = 0;  // x2 начальное

    // Определяем систему нелинейных уравнений (пересекающиеся окружности)
    DisjointCircles function; // Используется класс для описания окружностей

    // Создаем параметры системы
    SystemParameters params(n, m, maxIter, maxIterBeta, epsF, epsBeta, &function, x0);

    // Выбираем метод обработки (заменяем Symmetrization на Convolution или ExcludingRows)
    Convolution squaring; // Используется метод свертки

    // Создаем объект решения системы
    SystemOfNonlinearEquations solver(params, &squaring);

    // Пытаемся решить систему
    Vector solution = solver.Solve();

    // Выводим результат
    std::cout << "Решение системы:" << std::endl;
    std::cout << "x1 = " << solution(0) << std::endl;
    std::cout << "x2 = " << solution(1) << std::endl;

    // Анализируем результат (например, вывод нормы функции в конечной точке)
    double finalDiscrepancy = function.ComputeInPoint(solution).EuqlideanNorm();
    std::cout << "Норма функции в решении: " << finalDiscrepancy << std::endl;

    return 0;
}
