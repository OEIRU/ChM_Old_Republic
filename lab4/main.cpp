#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Типы данных
typedef vector<double> Vector;
typedef vector<vector<double>> Matrix;

// Функция для вычисления евклидовой нормы вектора
double euclideanNorm(const Vector &v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return sqrt(sum);
}

// Функция для умножения матрицы на вектор
Vector matrixVectorMultiply(const Matrix &mat, const Vector &vec) {
    Vector result(mat.size(), 0.0);
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[i].size(); ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}

// Функция для вычисления якобиана
Matrix formJacobiMatrix(const Vector &x, Vector F(const Vector &), int m, int n) {
    Matrix Jacobi(m, Vector(n, 0.0));
    double h = 1e-10;
    Vector temp = x;
    Vector Fp = F(x);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            temp[j] += h;
            Vector Fp_h = F(temp);
            for (int k = 0; k < m; ++k) {
                Fp_h[k] -= Fp[k];
            }
            Jacobi[i][j] = Fp_h[i] / h;
            temp[j] = x[j];
        }
    }
    return Jacobi;
}

// Функция для решения системы линейных уравнений методом Гаусса
Vector solveGauss(Matrix mat, Vector vec) {
    int n = vec.size();
    for (int i = 0; i < n; ++i) {
        double maxEl = abs(mat[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(mat[k][i]) > maxEl) {
                maxEl = abs(mat[k][i]);
                maxRow = k;
            }
        }

        swap(mat[maxRow], mat[i]);
        swap(vec[maxRow], vec[i]);

        for (int k = i + 1; k < n; ++k) {
            double c = -mat[k][i] / mat[i][i];
            for (int j = i; j < n; ++j) {
                if (i == j) {
                    mat[k][j] = 0;
                } else {
                    mat[k][j] += c * mat[i][j];
                }
            }
            vec[k] += c * vec[i];
        }
    }

    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = vec[i] / mat[i][i];
        for (int k = i - 1; k >= 0; --k) {
            vec[k] -= mat[k][i] * x[i];
        }
    }
    return x;
}

// Функция для свертки системы
void convolution(Matrix &mat, Vector &vec) {
    int m = mat.size();
    int n = mat[0].size();
    int rowsToDelete = m - n + 1;
    double sum = 0.0;
    Vector rowsSum(n, 0.0);

    for (int i = 0; i < rowsToDelete; ++i) {
        int row = 0;
        double min = fabs(vec[0]);
        for (int j = 0; j < m - i; ++j) {
            if (fabs(vec[j]) < min) {
                min = fabs(vec[j]);
                row = j;
            }
        }
        sum += pow(vec[row], 2);
        for (int j = 0; j < n; ++j) {
            rowsSum[j] += pow(mat[row][j], 2);
        }
        swap(vec[row], vec[m - i - 1]);
        swap(mat[row], mat[m - i - 1]);
    }
    vec[n - 1] = sum;
    for (int i = 0; i < n; ++i) {
        mat[n - 1][i] = rowsSum[i];
    }
}

// Функция для решения системы нелинейных уравнений методом Ньютона
Vector solveNewton(const Vector &x0, Vector F(const Vector &), int n, int m, int maxiter, double epsF) {
    Vector xk = x0;
    double F0Norm = euclideanNorm(F(x0));
    double FkNorm = F0Norm;

    for (int k = 0; k < maxiter && FkNorm / F0Norm > epsF; ++k) {
        Vector Fk = F(xk);
        Matrix Jacobi = formJacobiMatrix(xk, F, m, n);
        if (m != n) {
            convolution(Jacobi, Fk);
        }
        Vector dx = solveGauss(Jacobi, Fk);

        for (int i = 0; i < n; ++i) {
            xk[i] -= dx[i];
        }
        FkNorm = euclideanNorm(F(xk));

        cout << "Iteration " << k << ": x = ";
        for (double val : xk) {
            cout << val << " ";
        }
        cout << ", Discrepancy = " << FkNorm / F0Norm << endl;
    }
    return xk;
}

// Пример системы уравнений
Vector systemEquations(const Vector &x) {
    Vector F(2);
    F[0] = pow(x[0] + 2, 2) + pow(x[1] - 2, 2) - 4;
    F[1] = pow(x[0] - 2, 2) + pow(x[1] - 2, 2) - 4;
    return F;
}

// Основная функция
int main() {
    // Открываем файл для чтения входных данных
    ifstream inFile("input.txt");
    if (!inFile.is_open()) {
        cerr << "Не удалось открыть файл input.txt\n";
        return 1;
    }

    int n, m, maxiter;
    double epsF;
    Vector x0;

    // Считываем данные из файла
    inFile >> n >> m >> maxiter >> epsF;
    x0.resize(n);
    for (int i = 0; i < n; ++i) {
        inFile >> x0[i];
    }
    inFile.close();

    // Решение системы
    Vector solution = solveNewton(x0, systemEquations, n, m, maxiter, epsF);

    // Записываем результаты в файл
    ofstream outFile("output.txt");
    if (outFile.is_open()) {
        outFile << "Решение системы:\n";
        for (int i = 0; i < n; ++i) {
            outFile << "x" << i + 1 << " = " << fixed << setprecision(6) << solution[i] << "\n";
        }
        outFile.close();
    } else {
        cerr << "Не удалось открыть файл output.txt для записи.\n";
    }

    // Записываем итерации в файл для анализа
    ofstream iterFile("iterations.txt");
    if (iterFile.is_open()) {
        Vector xk = x0;
        double F0Norm = euclideanNorm(systemEquations(x0));
        double FkNorm = F0Norm;

        iterFile << "k,x1,x2,beta,discrepancy\n";
        for (int k = 0; k < maxiter && FkNorm / F0Norm > epsF; ++k) {
            Vector Fk = systemEquations(xk);
            Matrix Jacobi = formJacobiMatrix(xk, systemEquations, m, n);
            if (m != n) {
                convolution(Jacobi, Fk);
            }
            Vector dx = solveGauss(Jacobi, Fk);

            for (int i = 0; i < n; ++i) {
                xk[i] -= dx[i];
            }
            FkNorm = euclideanNorm(systemEquations(xk));

            iterFile << k << "," << xk[0] << "," << xk[1] << "," << 1.0 << "," << FkNorm / F0Norm << "\n";
        }
        iterFile.close();
    } else {
        cerr << "Не удалось открыть файл iterations.txt для записи.\n";
    }

    return 0;
}