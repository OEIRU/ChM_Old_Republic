#include "main.h"

template<typename T>
bool solveGauss(std::vector<std::vector<T>> A, std::vector<T>& x, std::vector<T> B)
{
    unsigned int n = A.size();
    // Проверка на согласованность размеров
    if (n == 0 || A[0].size() != n || B.size() != n) {
        return false; // Неверные размеры
    }

    for (unsigned int k = 0; k < n - 1; k++)
    {
        // Постолбцовый выбор главного элемента
        int m = k;
        for (int i = k + 1; i < n; i++)
        {
            if (std::abs(A[i][k]) > std::abs(A[m][k]))
            {
                m = i;
            }
        }
        if (std::abs(A[m][k]) < 1E-5)
        {
            return false; // Система вырождена
        }

        // Обмен строк
        if (m != k) {
            std::swap(A[k], A[m]);
            std::swap(B[k], B[m]);
        }

        // Прямой ход
        for (int i = k + 1; i < n; i++)
        {
            T factor = A[i][k] / A[k][k];
            B[i] -= factor * B[k];
            for (int j = k; j < n; j++)
            {
                A[i][j] -= factor * A[k][j];
            }
        }
    }

    // Обратный ход
    for (int k = n - 1; k >= 0; k--)
    {
        T sum = 0.0;
        for (int j = k + 1; j < n; j++)
        {
            sum += A[k][j] * x[j];
        }
        x[k] = (B[k] - sum) / A[k][k];
    }

    return true; // Успешное решение
}

