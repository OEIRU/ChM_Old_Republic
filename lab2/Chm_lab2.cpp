#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#define kol 14 
using namespace std;

typedef double type;

long n, m, n24, n15;
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

void fun1(type *x, type *xx) {
    long i2 = 0, i4 = 2, i5 = 1 + m + 2;
    updateValues(x, xx, 0, n15, i2, i4, i5, 1, 1 + m + 1);
    updateValues(x, xx, n15, n24, i2, i4, i5, 1, 1 + m + 1);
}

void fun2(type *x, type *xx) {
    long i1 = 0, i2 = 0, i4 = 2, i5 = 1 + m + 2; //allarma
    updateValues(x, xx, 0, i1, i2, i4, i5, 1, 1 + m + 1);
    updateValues(x, xx, i1, n15, i2, i4, i5, 1, 1 + m + 1);
    updateValues(x, xx, n15, n24, i2, i4, i5, 1, 1 + m + 1);
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
    from >> flag >> n >> m >> iter >> w >> e;
    n24 = n - 1;
    n15 = n - (1 + m + 1);

    di1 = new type[n15];
    di2 = new type[n24];
    di3 = new type[n];
    di4 = new type[n24];
    di5 = new type[n15];
    xx = new type[n];
    x = new type[n];
    f = new type[n];
    F = new type[n];

    for (long i = 0; i < n15; i++) from >> di1[i];
    for (long i = 0; i < n24; i++) from >> di2[i];
    for (long i = 0; i < n; i++) from >> di3[i];
    for (long i = 0; i < n24; i++) from >> di4[i];
    for (long i = 0; i < n15; i++) from >> di5[i];
    for (long i = 0; i < n; i++) from >> x[i];
    for (long i = 0; i < n; i++) from >> f[i];

    from.close();
}

type norma(type *A) {
    type s = 0;
    for (long i = 0; i < n; i++) {
        s += A[i] * A[i];
    }
    return sqrt(s);
}

void Ff() {
    fill(F, F + n, 0);
    long i1 = 1 + m + 1;
    long i2 = 1;

    for (int k = 0; k < n15; k++) {
        int I = k + i1;
        F[I] += di1[k] * x[k];
        F[k] += di5[k] * x[I];
    }
    for (int k = 0; k < n24; k++) {
        int I = k + i2;
        F[I] += di2[k] * x[k];
        F[k] += di4[k] * x[I];
    }
    for (int k = 0; k < n; k++) {
        F[k] += di3[k] * x[k];
    }
}

void Swap() {
    copy(xx, xx + n, x);
}

type otn_pogr_resh() {
    type *Xx = new type[n];
    for (int i = 0; i < n; i++) {
        Xx[i] = i + 1. - x[i];
    }
    type norm = norma(Xx);
    for (int i = 0; i < n; i++) {
        Xx[i] = 1. + i;
    }
    type Norm = norma(Xx);
    delete[] Xx;
    return norm / Norm;
}

type pogr(type NORM) {
    fill(F, F + n, 0);
    long i1 = 1 + m + 1;
    long i2 = 1;

    for (int k = 0; k < n15; k++) {
        int I = k + i1;
        F[I] += di1[k] * x[k];
        F[k] += di5[k] * x[I];
    }
    for (int k = 0; k < n24; k++) {
        int I = k + i2;
        F[I] += di2[k] * x[k];
        F[k] += di4[k] * x[I];
    }
    for (int k = 0; k < n; k++) {
        F[k] += di3[k] * x[k];
    }
    for (long i = 0; i < n; i++) {
        F[i] -= f[i];
    }
    return norma(F) / NORM;
}

int main() {
    clock_t t1, t2;
    t1 = clock();
    vvod();
    ofstream to("out.txt");
    ofstream to1("out1.txt");
    to.setf(ios_base::scientific, ios_base::floatfield);
    to.precision(kol);
    to1.setf(ios_base::scientific, ios_base::floatfield);
    to1.precision(2);

    type NORM = norma(f);
    type nevyazka = pogr(NORM);
    long I = 0;

    do {
        if (flag == 1) {
            (n15 <= 1 + m + 1) ? fun1(x, xx) : fun2(x, xx);
            Swap();
        } else {
            (n15 <= 1 + m + 1) ? fun1(x, x) : fun2(x, x);
        }

        vyvod(to, I++);
        nevyazka = pogr(NORM);
        type _x = otn_pogr_resh();
        to1 << "\n\ntek iter := " << I << "  otn.nevyazka := " << nevyazka << endl;
        to1 << "otn pogr resh := " << _x << endl;
        to1 << "число обусловленностей >=" << _x / nevyazka << "\n\n" << endl;

    } while (I < iter && nevyazka > e);

    t2 = clock();
    to1 << "время: " << (t2 - t1) / (double)CLOCKS_PER_SEC << " seconds";
    to.close();
    to1.close();
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
