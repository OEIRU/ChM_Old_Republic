#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

#define kol 14

using namespace std;

typedef double type;

type* di1;
type* di2;
type* di3;
type* di4;
type* di5;
type* f;
type* x;
type* xx;

type w, e;
int n, m, n24, n15, size_block, kol_block, iter, i1;

void LU_block() {
    int bcount = n / size_block;

    for (int k = 0; k < bcount; k++) {
        for (int i = 1; i < size_block; i++) {
            int l = (size_block - 1) * k + (i - 1) + k;
            di4[l] /= di3[(i - 1) + size_block * k];
            di3[i + size_block * k] -= di2[l] * di4[l];
        }
    }
}

type norma(const type* A) {
    type s = 0.;
    for (int i = 0; i < n; i++)
        s += A[i] * A[i];
    return sqrt(s);
}

void slau(bool F, int step) {
    int j = step * size_block;

    if (F) {
        for (int i = 0; i < size_block; i++)
            x[j + i] *= w;

        x[j] /= di3[j];

        for (int i = 1; i < size_block; i++) {
            x[j + i] -= di2[j + i - 1] * x[j + i - 1];
            x[j + i] /= di3[j + i];
        }

        for (int i = size_block - 2; i >= 0; i--)
            x[j + i] -= di4[j + i] * x[j + i + 1];
    } else {
        for (int i = 0; i < size_block; i++)
            x[j + i] += (1 - w) * xx[j + i];
    }
}

void out(ofstream& to, int I) {
    to.setf(ios_base::scientific, ios_base::floatfield);
    to.precision(kol);
    to << "\n\nIteration = " << I << ":\n\n";
    for (int i = 0; i < n; i++)
        to << x[i] << "\n";
    to << "Difference:\n\n";
    to.precision(2);
    for (int i = 0; i < n; i++)
        to << (i + 1. - x[i]) << "\n";
}

void enter() {
    ifstream from("in.txt");
    from >> n >> m >> size_block >> iter >> w >> e;

    kol_block = n / size_block;
    n24 = n - 1;
    n15 = n - (1 + m + 1);
    i1 = 1 + m + 1;

    di1 = new type[n15];
    di2 = new type[n24];
    di3 = new type[n];
    di4 = new type[n24];
    di5 = new type[n15];
    x = new type[n];
    xx = new type[n];
    f = new type[n];

    for (int i = 0; i < n15; i++) from >> di1[i];
    for (int i = 0; i < n24; i++) from >> di2[i];
    for (int i = 0; i < n; i++) from >> di3[i];
    for (int i = 0; i < n24; i++) from >> di4[i];
    for (int i = 0; i < n15; i++) from >> di5[i];
    for (int i = 0; i < n; i++) from >> x[i];
    for (int i = 0; i < n; i++) xx[i] = 0.;
    for (int i = 0; i < n; i++) from >> f[i];

    from.close();
}

void relax() {
    int nblock = i1 / size_block;
    int here = i1 % size_block;

    // Update x based on di5
    for (int i = 0; i < n15; i++)
        x[i] -= di5[i] * xx[i + i1];

    // Update last element of each block
    for (int i = 0; i < kol_block - 1; i++) {
        int j = (i + 1) * size_block - 1;
        x[j] -= di4[j] * xx[j + 1];
    }

    slau(true, 0);
    slau(false, 0);

    for (int i = 1; i < nblock; i++) {
        x[i * size_block] -= di2[i * size_block - 1] * x[i * size_block - 1];
        slau(true, i);
        slau(false, i);
    }

    int set = size_block - here;
    int j = nblock * size_block + here;

    for (int i = 0; i < set; i++)
        x[j + i] -= di1[i] * x[i];

    int i = nblock * size_block;
    x[i] -= di2[i - 1] * x[i - 1];

    slau(true, nblock);
    slau(false, nblock);

    for (int i = nblock + 1; i < kol_block; i++) {
        for (int j = 0; j < size_block; j++)
            x[i * size_block + j] -= di1[size_block * i - i1 + j] * x[size_block * i - i1 + j];

        x[i * size_block] -= di2[i * size_block - 1] * x[i * size_block - 1];

        slau(true, i);
        slau(false, i);
    }
}

void right_part(type* F, const type* X) {
    fill(F, F + n, 0.);
    for (int k = 0; k < n15; k++) {
        int I = k + i1;
        F[I] += di1[k] * X[k];
        F[k] += di5[k] * X[I];
    }
    for (int k = 0; k < n24; k++) {
        int I = k + 1;
        F[I] += di2[k] * X[k];
        F[k] += di4[k] * X[I];
    }
    for (int k = 0; k < n; k++)
        F[k] += di3[k] * X[k];
}

type norma(const type* x, const type* y) {
    type s = 0, s1 = 0;
    for (int i = 0; i < n; i++) {
        s += (x[i] - y[i]) * (x[i] - y[i]);
        s1 += y[i] * y[i];
    }
    return sqrt(s) / sqrt(s1);
}

type _pogr(const type* x, const type* F) {
    type* R = new type[n];
    right_part(R, x);
    type result = norma(R, F);
    delete[] R; // Free allocated memory
    return result;
}

int main() {
    clock_t tt1 = clock();
    ofstream to("out.txt");
    to.setf(ios_base::scientific, ios_base::floatfield);
    to.precision(kol);
    enter();

    type* F = new type[n];
    type* X = new type[n];
    for (int i = 0; i < n; i++)
        X[i] = i + 1.;

    LU_block();
    right_part(F, X);
    type t = 1;
    ofstream to1("out1.txt");
    to1.setf(ios_base::scientific, ios_base::floatfield);
    to1.precision(2);

    for (int I = 0; I < iter && t >= e; I++) {
        swap(xx, x); // Use std::swap if possible
        relax();
        out(to, I);

        t = norma(x, xx);
        type t1 = _pogr(x, F);
        to1 << "Iteration: " << I << endl;
        to1 << "Current residual: " << t << endl;
        to1 << "Current error: " << t1 << endl;
        to1 << "Condition number: " << t / t1 << "\n\n";
    }

    clock_t tt2 = clock();
    to1 << "Elapsed time: " << static_cast<double>(tt2 - tt1) / CLOCKS_PER_SEC << " seconds\n";
    to1.close();
    to.close();

    // Free allocated memory
    delete[] di1;
    delete[] di2;
    delete[] di3;
    delete[] di4;
    delete[] di5;
    delete[] f;
    delete[] x;
    delete[] xx;
    delete[] F;
    delete[] X;

    return 0;
}
