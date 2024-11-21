#include "main.h"

// Сгенерировать матрицу Гильберта
template<typename T>
void generateHilbert(SkylineStorageMatrix<T>& dst, std::vector<T>& dstr, int n)
{
    std::vector<std::vector<T>> square(n); // Для удобства параллельно профильной матрице
                                           // заполним квадратную
    for (int i = 0; i < n; i++)
    {
        square[i].resize(n);
    }

    // Матрица
    dst.ia.resize(n + 1);
    dst.ia[0] = 0;
    dst.ia[1] = 0;
    for (int i = 2; i < n + 1; i++)
    {
        dst.ia[i] = dst.ia[i - 1] + (i - 1);
    }

    int ii, jj;
    T sumij = (T)0;

    dst.al.resize(dst.ia[n]);
    ii = 1, jj = 0;
    for (int i = 0; i < dst.ia[n]; i++)
    {
        dst.al[i] = (T)((T)1 / (T)(ii + jj + 1));
        sumij += dst.al[i];

        square[ii][jj] = dst.al[i];
        if (++jj == ii)
        {
            ii++, jj = 0;
        }
    }

    dst.au.resize(dst.ia[n]);
    ii = 1, jj = 0;
    for (int i = 0; i < dst.ia[n]; i++)
    {
        dst.au[i] = (T)((T)1 / (T)(ii + jj + 1));
        sumij += dst.au[i];

        square[jj][ii] = dst.al[i];
        if (++jj == ii)
        {
            ii++, jj = 0;
        }
    }

    dst.di.resize(n);
    for (int i = 0; i < n; i++)
    {
        dst.di[i] = (T)((T)1 / (T)(i + i + 1));;
        square[i][i] = dst.di[i];
    }

    // Вектор
    dstr.resize(n);
    for (int i = 0; i < n; i++)
    {
        dstr[i] = (T)0;
        for (int j = 0; j < n; j++)
        {
            dstr[i] += square[i][j] * (T)(j + 1);
        }
    }

    // printMatrix(square); функция с записью в файл (пока не нужна)
    /*std::ostringstream os; os << n;
    std::fstream out("hilbert" + std::string(os.str()) + ".txt", std::ios::out);
    writeSkyline(dst, out);
    */


}

// Провести серию тестов с матрицами Гильберта размерности от n1 до n2
void hilbertTestSeries(int n1, int n2)
{
    // Открываем файл для записи
    std::ofstream outFile("hilbert_method.txt");
    if (!outFile.is_open()) {
        std::cerr << "Не удалось открыть файл для записи." << std::endl;
        return;
    }

    // Заголовок таблицы
    outFile << std::left << std::setw(5) << "k" 
            << std::setw(25) << "x^k (float)" 
            << std::setw(25) << "|x* - x^k| (float)" 
            << std::setw(25) << "x^k (double)" 
            << std::setw(25) << "|x* - x^k| (double)" << std::endl;
    outFile << std::string(105, '-') << std::endl;

    for (int n = n1; n <= n2; n++)
    {
        struct SkylineStorageMatrix<float> f_a;
        std::vector<float> f_b;
        generateHilbert(f_a, f_b, n);
        std::vector<float> f_x(n);
        bool ok1 = solve(f_x, f_a, f_b);

        struct SkylineStorageMatrix<double> d_a;
        std::vector<double> d_b;
        generateHilbert(d_a, d_b, n);
        std::vector<double> d_x(n);
        bool ok2 = solve(d_x, d_a, d_b);

        for (int i = 0; i < n; i++)
        {
                       outFile << std::left << std::setw(5) << (i == (n - 1) ? std::to_string(i + 1) : " ") 
                    << std::setw(25) << (ok1 ? std::to_string(f_x[i]) : "N/A")
                    << std::setw(25) << (ok1 ? std::to_string(fabs((i + 1) - f_x[i])) : "N/A")
                    << std::setw(25) << (ok2 ? std::to_string(d_x[i]) : "N/A")
                    << std::setw(25) << (ok2 ? std::to_string(fabs((i + 1) - d_x[i])) : "N/A")
                    << std::endl;
        }

        // Разделительная линия между таблицами для разных n
        if (n != n2) {
            outFile << std::string(105, '-') << std::endl;
        }
    }

    // Закрываем файл
    outFile.close();
}
