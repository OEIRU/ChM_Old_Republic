#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace std::chrono;

static void readInf(const string& filename, int& n, int& m, double& omega, int& MAXITER, double& tol, int& flag);
static double sum(const vector<vector<double>>& A, const vector<double>& y, const int& pos, const int& m, const int& flag);
static double getNorm(const vector<double>& F);
static double getFinalNorm(const vector<vector<double>>& A, const vector<double>& F, const vector<double>& x, const int& m);
static void Yakoby_Zaydelya(const vector<vector<double>>& A, const vector<double>& b, vector<double>& X, vector<double>& Xprev, const int& m, const double& omega,
   int& MAXITER, double& tol, const bool& flag);
static void readVect(const string& filename, vector<double>& vec, const int& n);
static void readMatrix(string& filename, vector<vector<double>>& A, const int& n);
static void writeVectorToFile(vector<double>& x, const string& filename);


int main() {
   setlocale(LC_ALL, "RUS");
   string file_matrix = "matrix.txt";
   string file_vector = "vector.txt";
   string file_index = "index.txt";
   string file_x_vect = "x_vector.txt";
   string file_x_new_vect = "output.txt";
   int n, m, flag, MAXITER;
   double omega, tol;
   vector<vector<double>> A;
   vector<double> F;
   vector<double> x;
   vector<double> x_prev;
   vector<double> x_toch;
   vector<double> x_dif;
   readInf(file_index, n, m, omega, MAXITER, tol, flag);
   x_toch.resize(n);
   x_dif.resize(n);
   F.resize(n);
   x.resize(n);
   A.resize(5);
   for (auto& row : A) row.resize(n);
   readMatrix(file_matrix, A, n);
   readVect(file_vector, F, n);
   readVect(file_x_vect, x, n);
   for (int i = 0; i < n; i++) {
      x_toch[i] = i + 1;
   }
   if (flag == 0) {
      Yakoby_Zaydelya(A, F, x, x_prev, m, omega, MAXITER, tol, flag);
   }
   else {
      Yakoby_Zaydelya(A, F, x, x, m, omega, MAXITER, tol, flag);
   }

   /*for (double w = 0; w <= 1.1; ) {
      cout << "w = " <<  w << endl;
      Yakoby_Zaydelya(A, F, x, x, m, w, MAXITER, tol, flag);
      for (int i = 0; i < n; i++) {
         x_dif[i] = x_toch[i] - x[i];
         cout << x_toch[i] - x[i] << endl;
                                                                                       // код для исследований
      }
      cout << "--------" << endl;
      cout << (getNorm(x_dif) / getNorm(x_toch)) / getFinalNorm(A, F, x, m) << endl;
      cout << "--------" << endl;
      for (int i = 0; i < n; i++) {
         cout << x[i] << endl;
      }
      readVect(file_x_vect, x, n);
      w += 0.1;
      cout << "--------" << endl;
   }*/
   writeVectorToFile(x, file_x_new_vect);

   return 0;
}

static void readInf(const string& filename, int& n, int& m, double& omega, int& MAXITER, double& tol, int& flag)
{
   ifstream fin(filename);

   if (!fin)
   {
      cerr << "Error reading pars file." << endl;
      exit(1);
   }

   fin >> n >> m >>omega >> MAXITER >> tol >> flag;
   cout << n << " " << m << " " << omega << " " << MAXITER << " " << tol << " " << flag<<endl;
   fin.close();
}
static void readMatrix(string& filename, vector<vector<double>>& A, const int& n) {
   ifstream fin(filename);
   if (!fin)
   {
      cerr << "Error reading matrix file." << endl;
      exit(1);
   }
   for (int i = 0; i != 5; i++) {
      for (int j = 0; j != n; j++) {
         fin >> A[i][j];
      }
   }

   fin.close();
}
static void readVect(const string& filename, vector<double>& vec, const int& n)
{
   ifstream fin(filename);

   if (!fin)
   {
      cerr << "Error reading vector file." << endl;
      exit(1);
   }

   for (int i = 0; i < n; i++)
      fin >> vec[i];

   fin.close();
}
static double getNorm(const vector<double>& F)
{
   int n = int(F.size());

   double summ = 0;

   if (n < 1) return 0;

   for (int i = 0; i != n; i++)
      summ += F[i] * F[i];

   return sqrt(summ);
}
static double getFinalNorm(const vector<vector<double>>& A, const vector<double>& F, const vector<double>& x, const int& m)
{
   vector<double> vec;
   int n = (int)x.size();
   vec.resize(n);

   for (int i = 0; i < n; i++)
      vec[i] = F[i] - sum(A, x, i, m, 1);

   return getNorm(vec) / getNorm(F);
}
static double sum(const vector<vector<double>>& A, const vector<double>& y, const int& pos, const int& m, const int& flag)
{
   double summ = 0;
   if (pos > 0) summ += A[1][pos] * y[pos - 1];
   if (flag == 1 || flag == 0) summ += A[2][pos] * y[pos];
   if (pos > m + 1) summ += A[0][pos] * y[pos - m - 2];
   if (pos < y.size() - 1) summ += A[3][pos] * y[pos + 1];
   if (pos < y.size() - 2 - m) summ += A[4][pos] * y[pos + m + 2];

   return summ;
}
static void Yakoby_Zaydelya(const vector<vector<double>>& A, const vector<double>& b, vector<double>& X, vector<double>& Xprev, const int& m, const double& omega,
   int& MAXITER, double& tol, const bool& flag)
{
   int n = (int)b.size();
   double summ;
   double norm = 1;
   int i;
   for (i = 0; i != MAXITER && norm >= tol; i++)
   {
      Xprev = X;
      for (int j = 0; j < n; j++)
      {
         summ = sum(A, Xprev, j, m, flag);
         X[j] = Xprev[j] + omega / A[2][j] * (b[j] - summ);
      }
      norm = getFinalNorm(A, b, X, m);
      //cout << "iteration number = " << i << endl; // Условие
      //cout << "norm = " << norm << endl; // Условие
   }

   if (norm >= tol)
      cout << "exit by max_iterations" << endl;
   else
      cout << "exit by " << i << " iterations" << endl;

}
static void writeVectorToFile(vector<double>& x, const string& filename) {
   ofstream file(filename);

   if (!file.is_open()) {
      cerr << "Error reading output file." << endl;
      return;
   }

   for (double value : x) {
      file << scientific << setprecision(15) << value << "\n";
   }

   file.close();
}