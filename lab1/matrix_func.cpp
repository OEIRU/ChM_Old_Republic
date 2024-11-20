#include "main.h"

template<typename T>
void writeVector(std::vector<T>& from, std::ostream& to)
{
    to << from.size() << std::endl;

    for (unsigned int i = 0; i < from.size(); i++)
    {
        to << std::fixed << std::setprecision(14) << from[i] << std::endl;
    }
    to << std::endl;
}

template<typename T>
void writeSkyline(SkylineStorageMatrix<T>& from, std::ostream& to)
{
    to << from.di.size() << std::endl;
    for (unsigned int i = 0; i < from.ia.size(); i++)
    {
        to << from.ia[i] << " ";
    }
    to << std::endl;

    for (unsigned int i = 0; i < from.al.size(); i++)
    {
        to << std::fixed << std::setprecision(14) << from.al[i] << " ";
    }
    to << std::endl;
    for (unsigned int i = 0; i < from.au.size(); i++)
    {
        to << std::fixed << std::setprecision(14) << from.au[i] << " ";
    }
    to << std::endl;

    for (unsigned int i = 0; i < from.di.size(); i++)
    {
        to << std::fixed << std::setprecision(14) << from.di[i] << " ";
    }
    to << std::endl;
}


 // Считать вектор
template<typename T>
void readVector(std::istream& from, std::vector<T>& to)
{
    int n;
    from >> n;
    to.resize(n);
    for (int i = 0; i < n; i++)
    {
        from >> to[i];
    }
}

// Распечатать квадратную матрицу
template<typename T>
void printMatrix(const std::vector<std::vector<T>> a)
{
    for (unsigned int i = 0; i < a.size(); i++)
    {
        for (unsigned int j = 0; j < a.size(); j++)
        {
            std::cout << std::setw(20) << a[i][j];
        }
        std::cout << std::endl;
    }
}

// Распечатать матрицу в профильном формате
template<typename T>
void printSkyline(SkylineStorageMatrix<T>& a)
{
    std::vector<std::vector<T>> s;
    skylineToSquare(a, s);
    printMatrix(s);
}


// LU(sq)-разложение матрицы в профильном формате
template<typename T, typename T2>
bool LusqSkyline(SkylineStorageMatrix<T>& a)
{
    // Первая строка матрицы
    _sqrt(a.di[0], a.di[0]); // d_1 = sqrt(a_11);
    for (unsigned int i = 1; i < a.di.size(); i++)
    {
        // Нужны ненулевые элементы первой строки
        // i'ый элемент первой строки ненулевой <=> i'ый столбец содержит i элементов
        if (elementsInLine(i) == i)
        {
            a.au[a.ia[i]] = a.au[a.ia[i]] / a.di[0]; // U_1i = a_1i/q1
        }
    }

    // Обход по строкам матрицы
    for (unsigned int i = 1; i < a.di.size(); i++)
    {
        if (a.al.size())
        {
            // l_ij = (a_ij - \sum_{k=1}{j-1} l_ik * u_kj) / d_j
            /* Иными словами, для вычисления j'го элемента i'ой строки, необходимо найти
            * сумму произведений первых j-1 элементов i'ой строки на
            * соответствующие первые j-1 элементов j'го столбца.
            * Затем вычесть её из a_ij и полученную разность поделить на a_jj. Легко
            * заметить, что подобное преобразование сохраняет портрет матрицы.
            */
            unsigned int alStartIndex = a.ia[i]; // Индекс первого ненулевого элемента i'ой строки в массиве al
            unsigned int alEndIndex = a.ia[i] + elementsInLine(i); // До сюда идти в цикле
            unsigned int j = i - elementsInLine(i);
            // До первого ненулевого элемента не было ненулевых элементов
            a.al[alStartIndex] = a.al[alStartIndex] / a.di[j];
            alStartIndex++;
            j++;
            // Основной цикл
            for (unsigned int alIndex = alStartIndex; alIndex < alEndIndex; alIndex++, j++)
            {
                T sum = (T)0;
                unsigned int elementsInLineI = elementsInLine(i);
                unsigned int elementsInColJ = elementsInLine(j);
                for (unsigned int k = 0; k < j; k++)
                {
                    T l_ik = (k < i - elementsInLineI) ? (T)0 : a.al[a.ia[i] + (k - (i - elementsInLineI))];
                    T u_kj = (k < j - elementsInColJ) ? (T)0 : a.au[a.ia[j] + (k - (j - elementsInColJ))];
                    sum += l_ik * u_kj;
                }
                a.al[alIndex] = (a.al[alIndex] - sum) / a.di[j];
            }
        }

        // d_i = sqrt(d_i - \sum_{k=1}{i-1} l_ik * u_ki)
        {
            T2 sum = (T2)0;
            unsigned int elementsInLineI = elementsInLine(i);
            for (unsigned int k = 0; k < i; k++)
            {
                T l_ik = (k < i - elementsInLineI) ? (T)0 : a.al[a.ia[i] + (k - (i - elementsInLineI))];
                T l_ki = (k < i - elementsInLineI) ? (T)0 : a.au[a.ia[i] + (k - (i - elementsInLineI))];
                sum += (T2)l_ik * (T2)l_ki;
            }
            // printf("%d %lf\n", i, sum);
            _sqrt(a.di[i], a.di[i] - sum);
        }

        if (a.au.size())
        {
            // u_ij = (a_ij - \sum_{k=1}{i-1} l_ik * u_kj) / d_j
            /* То же самое, что l_ij, только писать нужно в al (причём, не последовательно) и сумма по k до i,
            * и делить на di[i]
            */
            for (unsigned int j = i + 1; j < a.di.size(); j++)
            {
                if (i >= j - elementsInLine(j))
                {
                    T2 sum = (T2)0;
                    unsigned int elementsInLineI = elementsInLine(i);
                    unsigned int elementsInColJ = elementsInLine(j);
                    for (unsigned int k = 0; k < i; k++)
                    {
                        T l_ik = (k < i - elementsInLineI) ? (T)0 : a.al[a.ia[i] + (k - (i - elementsInLineI))];
                        T u_kj = (k < j - elementsInColJ) ? (T)0 : a.au[a.ia[j] + (k - (j - elementsInColJ))];
                        sum += (T2)l_ik * (T2)u_kj;
                    }

                    unsigned int auIndex = a.ia[j] + (i - (j - elementsInLine(j)));
                    a.au[auIndex] = (a.au[auIndex] - sum) / a.di[i];
                }
            }
        }
    }
    // printSkyline(a);
    return true;
}

