#include "main.h"

// Сравнить Гаусса и разложение тестами на числом обусловленности
void compareGaussLusq(int tests)
{
    int counter = 1; 
    std::ofstream outFile("сomparison_methods.txt");
    if (!outFile.is_open()) {
        std::cerr << "Не удалось открыть файл для записи." << std::endl;
        return;
    }

    // Заголовок таблицы
    outFile << std::left
              << std::setw(5)  << "k"
              << std::setw(25) << "x^k (LU(sq))"
              << std::setw(20) << "|x* - x^k| (LU(sq))"
              << std::setw(25) << "x^k (Gauss)"
              << std::setw(20) << "|x* - x^k| (Gauss)"
              << std::endl;
    outFile << std::string(120, '-') << std::endl;

    for (int test = 0; test < tests; test++)
    {
        std::fstream d_ifile("test6.txt");
        struct SkylineStorageMatrix<double> d_a;
        readSkyline(d_ifile, d_a);
        int n = d_a.di.size();
        std::vector<double> d_b;
        readVector(d_ifile, d_b);
        d_a.di[0] += (double)pow(10.0, -test);
        d_b[0] += (double)(pow(10.0, -test));
        std::vector<double> d_x(n);
        bool ok1 = solve(d_x, d_a, d_b);

        std::fstream g_ifile("test6.txt");
        struct SkylineStorageMatrix<double> g_a;
        readSkyline(g_ifile, g_a);
        std::vector<double> g_b;
        readVector(g_ifile, g_b);
        g_a.di[0] += (double)pow(10.0, -test);
        g_b[0] += (double)(pow(10.0, -test));
        std::vector<std::vector<double>> A;
        skylineToSquare(g_a, A);
        std::vector<double> g_x(n);
        bool ok2 = solveGauss(A, g_x, g_b);

        // Печать строк таблицы
        for (int i = 0; i < n; i++)
        {
                outFile << std::left << std::setw(5) << (counter%10 == 0? std::to_string(counter/10) : " ")
            << std::setw(25) << (ok1 ? std::to_string(d_x[i]) : "N/A")
            << std::setw(20) << (ok1 ? std::to_string(fabs((i + 1) - d_x[i])) : "N/A")
            << std::setw(25) << (ok2 ? std::to_string(g_x[i]) : "N/A")
            << std::setw(20) << (ok2 ? std::to_string(fabs((i + 1) - g_x[i])) : "N/A")
            << std::endl;
            counter++;
        }

        // Разделительная линия между тестами
        if (test != tests - 1) {
            outFile << std::string(120, '-') << std::endl;
        }
    }
}

// Провести серию тестов на число обусловленности
void conditionNumberTestSeries(int tests)
{
    int counter = 1;
    // Открываем файл для записи
    std::ofstream outFile("condition_method.txt");
    if (!outFile.is_open()) {
        std::cerr << "Не удалось открыть файл для записи." << std::endl;
        return;
    }

    // Заголовок таблицы
    outFile << std::left 
              << std::setw(5)  << "k" 
              << std::setw(25) << "x^k (float)" 
              << std::setw(20) << "|x* - x^k| (float)" 
              << std::setw(25) << "x^k (double)" 
              << std::setw(20) << "|x* - x^k| (double)" 
              << std::setw(25) << "x^k (mixed)" 
              << std::setw(20) << "|x* - x^k| (mixed)" 
              << std::endl;
    outFile << std::string(140, '-') << std::endl;

    for (int test = 0; test < tests; test++)
    {
        // Подготовка данных для float
        std::fstream f_ifile("test6.txt");
        struct SkylineStorageMatrix<float> f_a;
        readSkyline(f_ifile, f_a);
        int n = f_a.di.size();
        std::vector<float> f_b;
        readVector(f_ifile, f_b);
        f_a.di[0] += (float)pow(10.0, -test);
        f_b[0] += (float)pow(10.0, -test);
        std::vector<float> f_x(n);
        bool ok1 = solve(f_x, f_a, f_b);

        // Подготовка данных для double
        std::fstream d_ifile("test6.txt");
        struct SkylineStorageMatrix<double> d_a;
        readSkyline(d_ifile, d_a);
        std::vector<double> d_b;
        readVector(d_ifile, d_b);
        d_a.di[0] += (double)pow(10.0, -test);
        d_b[0] += (double)pow(10.0, -test);
        std::vector<double> d_x(n);
        bool ok2 = solve(d_x, d_a, d_b);

        // Подготовка данных для mixed
        std::fstream m_ifile("test6.txt");
        struct SkylineStorageMatrix<float> m_a;
        readSkyline(m_ifile, m_a);
        std::vector<float> m_b;
        readVector(m_ifile, m_b);
        m_a.di[0] += (float)pow(10.0, -test);
        m_b[0] += (float)pow(10.0, -test);
        std::vector<float> m_x(n);
        bool ok3 = solve<float, double>(m_x, m_a, m_b);

        // Печать строк таблицы
        for (int i = 0; i < n; i++)
        {
            outFile << std::left << std::setw(5) << (counter%10 == 0? std::to_string(counter/10) : " ")
                      << std::setw(25) << (ok1 ? std::to_string(f_x[i]) : "N/A")
                      << std::setw(20) << (ok1 ? std::to_string(fabs((i + 1) - f_x[i])) : "N/A")
                      << std::setw(25) << (ok2 ? std::to_string(d_x[i]) : "N/A")
                      << std::setw(20) << (ok2 ? std::to_string(fabs((i + 1) - d_x[i])) : "N/A")
                      << std::setw(25) << (ok3 ? std::to_string(m_x[i]) : "N/A")
                      << std::setw(20) << (ok3 ? std::to_string(fabs((i + 1) - m_x[i])) : "N/A")
                      << std::endl;
            counter++;
        }

        // Разделительная линия между тестами
        if (test != tests - 1) {
            outFile << std::string(140, '-') << std::endl;
        }
    }
    outFile.close();
}