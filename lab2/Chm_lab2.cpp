#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#define kol 14 
using namespace std;

typedef double type;

long n, m;
int iter, flag;
type w, e;
type *x, *xx, *f, *F;
type *di1, *di2, *di3, *di4, *di5;

void updateValues(type *x, type *xx, long start, long end, long &i2, long &i4, long &i5, long j4, long j5) {
    for (long i = start; i < end; i++, i2++, i4++, i5++) {
        type s = f[i];
        s -= di3[i] * x[i] + di4[i] * x[i4] + di5[i] * x[i5];
        s = w * s / di3[i];
        xx[i] = x[i] + s;
    }
}

void Jacobi(type *x, type *xx) {
    long i2 = 0, i4 = 2, i5 = 1 + m + 2;
    updateValues(x, xx, 0, n, i2, i4, i5, 1, 1 + m + 1);
}

void GaussSeidel(type *x, type *xx) {
    long i2 = 0, i4 = 2, i5 = 1 + m + 2;
    updateValues(x, xx, 0, n, i2, i4, i5, 1, 1 + m + 1);
    for (long i = 0; i < n; i++) {
        x[i] = xx[i];
    }
}

void BlockRelaxation(type *x, type *xx) {
    // Implement the block relaxation method here
}

void vyvod(ofstream &to, int i) {
    to.setf(ios_base::scientific, ios_base::floatfield);
    to.precision(kol);
    to << "\n\n итерация = " << i << ":\n\n";
    for (long i = 0; i < n; i++) {
        to << x[i] << "\n";
    }
    to << "разность:\n\n";
    to.setf(ios_base::scientific, ios_base::floatfield);
    to.precision(2);
    for (long i = 0; i < n; i++) {
        to << i + 1. - x[i] << "\n";
    }
}

void vvod() {
    ifstream from("in.txt");
    if (!from) {
        cerr << "Ошибка открытия файла!" << endl;
        exit(1);
    }

    from >> flag >> n >> m >> iter >> w >> e;

    // Initialize arrays
    di1 = new type[n];
    di2 = new type[n];
    di3 = new type[n];
    di4 = new type[n];
    di5 = new type[n];
    xx = new type[n];
    x = new type[n];
    f = new type[n];
    F = new type[n];

    // Read data from file
    for (long i = 0; i < n; i++) from >> di1[i]; // Read first array
    for (long i = 0; i < n; i++) from >> di2[i]; // Read second array
    for (long i = 0; i < n; i++) from >> di3[i]; // Read third array
    for (long i = 0; i < n; i++) from >> di4[i]; // Read fourth array
    for (long i = 0; i < n; i++) from >> di5[i]; // Read fifth array
    for (long i = 0; i < n; i++) from >> x[i]; // Initial approximation
    for (long i = 0; i < n; i++) from >> f[i]; // Right-hand side vector

    from.close();
}

type norma(type *A) {
    type s = 0;
    for (long i = 0; i < n; i++) {
        s += A[i] * A[i];
    }
    return sqrt(s);
}

void calculateResidual() {
    fill(F, F + n, 0);
    for (int k = 0; k < n; k++) {
        F[k] += di1[k] * x[k];
        F[k] += di2[k] * x[k];
        F[k] += di3[k] * x[k];
        F[k] += di4[k] * x[k];
        F[k] += di5[k] * x[k];
    }
    for (long i = 0; i < n; i++) {
        F[i] -= f[i];
    }
}

type calculateRelativeResidual() {
    calculateResidual();
    return norma(F);
}

int main() {
    clock_t t1, t2;
    t1 = clock();
    vvod();
    ofstream to("out.txt");
    to.setf(ios_base::scientific, ios_base::floatfield);
    to.precision(kol);

    type nevyazka;
    long I = 0;

    do {
        if (flag == 1) {
            Jacobi(x, xx);
        } else if (flag == 2) {
            GaussSeidel(x, xx);
        } else if (flag == 3) {
            BlockRelaxation(x, xx);
        }

        vyvod(to, I++);
        nevyazka = calculateRelativeResidual();
        cout << "Текущая итерация: " << I << ", относительная невязка: " << nevyazka << endl;

    } while (I < iter && nevyazka > e);

    t2 = clock();
    to << "время: " << (t2 - t1) / (double)CLOCKS_PER_SEC << " seconds";
    to.close();
    delete[] di1;
    delete[] di2;
    delete[] di3;
    delete[] di4;
    delete[] di5;
    delete[] xx;
    delete[] x;
    delete[] f;
    delete[] F;

    return 0;
}
