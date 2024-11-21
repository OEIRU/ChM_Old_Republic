#include "main.h"

// Считать матрицу в профильном формате
template<typename T>
void readSkyline(std::istream& from, SkylineStorageMatrix<T>& to)
{
    int n;
    from >> n;
    to.ia.resize(n + 1);
    for (int i = 0; i < n + 1; i++)
    {
        from >> to.ia[i];
    }

    int k = to.ia[n];
    to.al.resize(k);
    for (int i = 0; i < k; i++)
    {
        from >> to.al[i];
    }
    to.au.resize(k);
    for (int i = 0; i < k; i++)
    {
        from >> to.au[i];
    }

    to.di.resize(n);
    for (int i = 0; i < n; i++)
    {
        from >> to.di[i];
    }
}

// Сконвертировать матрицу в профильном формате в квадратную
template<typename T>
void skylineToSquare(SkylineStorageMatrix<T>& a, std::vector<std::vector<T>>& s)
{
    #define elementsInLine(i) (a.ia[i + 1] - a.ia[i])

    s.resize(a.di.size());
    for (unsigned int i = 0; i < a.di.size(); i++)
    {
        s[i].resize(a.di.size());
        for (unsigned int j = 0; j < a.di.size(); j++)
        {
            s[i][j] = (T)0;
        }
    }

    for (unsigned int i = 0; i < a.di.size(); i++)
    {
        s[i][i] = a.di[i];
    }

    for (unsigned int i = 1; i < a.di.size(); i++)
    {
        int offset = i - elementsInLine(i);
        for (int j = 0; j < elementsInLine(i); j++)
        {
            s[i][offset + j] = a.al[a.ia[i] + j];
            s[offset + j][i] = a.au[a.ia[i] + j];
        }
    }

    #undef elementsInLine
}

