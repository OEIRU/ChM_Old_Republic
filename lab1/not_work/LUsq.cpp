#include "main.h"


// Решить уравнение Ly=b
template<typename T, typename T2>
bool findY(std::vector<T>& y, SkylineStorageMatrix<T>& a, std::vector<T>& b)
{
    /* Для удобства. Присвоить dst значение частного a и b.
     * Если b = 0, вернуть false из текущей подпрограммы.
     */
    //#define _div(dst, a, b) { if (((b) > (T)0 ? (b) : -(b)) < (T)EPS) return false; dst = (a) / (b); }

    // Количество элементов в i'ой строке/столбце
    //#define elementsInLine(i) (a.ia[i + 1] - a.ia[i])

    for (unsigned int i = 0; i < a.di.size(); i++)
    {
        T2 sum = (T2)0;
        unsigned int elementsInLineI = elementsInLine(i);
		for (unsigned int alIndex = a.ia[i], bIndex = i - elementsInLineI, counter = 0; counter < elementsInLineI; alIndex++, bIndex++, counter++)
		{
			sum += (T2)a.al[alIndex] * (T2)y[bIndex];
		}

        _div(y[i], b[i] - sum, a.di[i]);
    }

    //#undef _div
	//#undef elementsInLine
    return true;
}

// Решить уравнение Ux=y
template<typename T, typename T2>
bool findX(std::vector<T>& x, SkylineStorageMatrix<T>& a, std::vector<T>& y)
{
    /* Для удобства. Присвоить dst значение частного a и b.
     * Если b = 0, вернуть false из текущей подпрограммы.
     */
    //#define _div(dst, a, b) { if (((b) > (T)0 ? (b) : -(b)) < (T)EPS) return false; dst = (a) / (b); }

    // Количество элементов в i'ой строке/столбце
    //#define elementsInLine(i) (a.ia[i + 1] - a.ia[i])

	for (unsigned int i = a.di.size() - 1; true; i--)
	{
        T2 sum = (T2)0;
        for (unsigned int j = i + 1; j < a.di.size(); j++)
		{
			unsigned int elementsInColJ = elementsInLine(j);
            sum += (i < j - elementsInColJ) ? (T2)0 : (T2)a.au[a.ia[j] + (i - (j - elementsInColJ))] * x[j];
		}

        _div(x[i], y[i] - sum, a.di[i]);
		if (i == 0) break;
	}

    //#undef _div
	//#undef elementsInLine
    return true;
}

// Решить уравнение Ax=B
template<typename T, typename T2>
bool solve(std::vector<T>& x, SkylineStorageMatrix<T>& a, std::vector<T>& b)
{
    //#define _div(dst, a, b) { if (((b) > (T)0 ? (b) : -(b)) < (T)EPS) return false; dst = (a) / (b); }
    //#define elementsInLine(i) (a.ia[i + 1] - a.ia[i])

	std::vector<T> y(b.size());

    if (!LusqSkyline<T,T2>(a))
    {
        std::cerr << "Matrix is not LU(sq)-decomposeable" << std::endl;
        return false;
    }

    if (!findY<T,T2>(y, a, b) || !findX<T,T2>(x, a, y))
    {
        std::cerr << "System is inconsistent" << std::endl;
        return false;
    }
    //#undef _div
	//#undef elementsInLine

	return true;
}
