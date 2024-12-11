#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <math.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

void Write_toFile(vector<double>x)
{
   ofstream outFile("output.txt");

   if (!outFile.is_open()) {
      cerr << "Не удалось открыть файл для записи!" << endl;
   }

   for (const auto& element : x) {
      outFile << scientific << setprecision(15) << element << endl;
   }

   outFile.close();
}

void readInf(const string& filename, int& N, int& MAXITER, double& tol)
{
   ifstream fin(filename);

   if (!fin)
   {
      cerr << "Error reading pars file." << endl;
      exit(1);
   }

   fin >> N >> MAXITER >> tol;
   fin.close();
}
void Initilize(vector<int>& ig, vector<double>& di, vector<double>& ggl, vector<int>& jg, vector<double>& pr, const int& n)
{
   ifstream Ggl("GGL.txt");
   ifstream Jg("JG.txt");
   ifstream Ig("IG.txt");
   ifstream Pr("PR.txt");
   ifstream Di("DI.txt");
   ig.resize(n + 1);
   for (int i = 0; i < n + 1; i++)
   {
      Ig >> ig[i];
   }

   di.resize(n);
   for (int i = 0; i < n; i++)
   {
      Di >> di[i];
   }

   ggl.resize(ig[n] - 1);
   for (int i = 0; i < ig[n] - 1; i++)
   {
      Ggl >> ggl[i];
   }
   jg.resize(ig[n] - 1);
   for (int i = 0; i < ig[n] - 1; i++)
   {
      Jg >> jg[i];
   }
   pr.resize(n);
   for (int i = 0; i < n; i++)
   {
      Pr >> pr[i];
   }

}
double ScalMult(vector<double> x, vector<double> y)
{
   double summ = 0;
   for (int i = 0; i < x.size(); i++) {
      summ += x[i] * y[i];
   }
   return summ;
}
vector<double> Mult(vector<double>& di, vector<double>& ggl, vector<int>& jg, vector<int>& ig, int& N, vector<double>& x)
{
   vector<double> y;
   y.resize(N);
   int j = 0;
   for (int i = 0; i < N; i++) {
      y[i] += di[i] * x[i];
      for (int k = 0; k < ig[i + 1] - ig[i]; k++) {
         y[i] += ggl[j] * x[jg[j]];
         y[jg[j]] += ggl[j] * x[i];
         j++;
      }
   }

   return y;
}
vector<double> Mult(vector<double>& di, int& N, vector<double>& x)
{
   vector<double> y;
   y.resize(N);
   int j = 0;
   for (int i = 0; i < N; i++) {
      y[i] += sqrt((1 / di[i])) * x[i];
   }
   return y;
}

double getNorm(const vector<double>& x)
{
   int n = int(x.size());

   double summ = 0;

   if (n < 1) return 0;

   for (int i = 0; i != n; i++)
      summ += x[i] * x[i];

   return sqrt(summ);
}

void LOS(int& N, vector<double>& di, vector<double>& ggl, vector<int>& jg, vector<double>& pr, vector<int>& ig, int& MAXITER, double& tol)
{
   double norm_square = 1.0;
   int j = 0;
   int i;
   double alpha;
   double beta;
   double r_old;
   double r_too_old;
   double p_old;
   vector<double> r;
   vector<double> x;
   vector<double> z;
   vector<double> p;
   vector <double> Ap;
   r.resize(N);
   x.resize(N);
   z.resize(N);
   p.resize(N);
   for (i = 0; i < N; i++) {
      r[i] += di[i] * x[i];
      for (int k = 0; k < ig[i + 1] - ig[i]; k++) {
         r[i] += ggl[j] * x[jg[j]];
         r[jg[j]] += ggl[j] * x[i];
         j++;
      }
      r[i] = pr[i] - r[i];
   }
   z = r;
   p = Mult(di, ggl, jg, ig, N, z);
   for (i = 0; i != MAXITER && norm_square >= tol; i++)
   {
      r_too_old = ScalMult(r, r);
      r_old = ScalMult(p, r);
      p_old = ScalMult(p, p);
      alpha = r_old / p_old;
      for (j = 0; j < x.size(); j++)
      {
         x[j] = x[j] + alpha * z[j];
         r[j] = r[j] - alpha * p[j];
      }
      Ap = Mult(di, ggl, jg, ig, N, r);
      beta = -1 * ScalMult(p, Ap) / p_old;
      for (j = 0; j < x.size(); j++)
      {
         z[j] = r[j] + beta * z[j];
         p[j] = Ap[j] + beta * p[j];
      }
      for (j = 0; j < r.size(); j++) {
         norm_square = r_too_old - pow(alpha, 2) * p_old; // Использовать предыдущую невязку
      }
      cout << "iterations = " << i+1<< endl;
      cout << "norm_square = " << norm_square << endl;
   }
   if (norm_square >= tol) cout << "exit by max_iterations" << endl;
   else cout << "exit by " << i << " iterations" << endl;
   /*for (i = 0; i < x.size(); ++i) {
      cout << scientific << setprecision(15) << x[i] << " ";
   }*/
   Write_toFile(x);
}

void LOS_diag(int& N, vector<double>& di, vector<double>& ggl, vector<int>& jg, vector<double>& pr, vector<int>& ig, int& MAXITER, double& tol)
{
   double norm_square = 1.0;
   int j = 0;
   int i;
   double alpha;
   double beta;
   double pr_old;
   double r_old;
   double p_old;
   vector<double> temp;
   vector<double> temp1;
   vector<double> r;
   vector<double> x;
   vector<double> z;
   vector<double> p;
   temp.resize(N);
   r.resize(N);
   x.resize(N);
   z.resize(N);
   p.resize(N);
   for (i = 0; i < N; i++) {
      r[i] += di[i] * x[i];
      for (int k = 0; k < ig[i + 1] - ig[i]; k++) {
         r[i] += ggl[j] * x[jg[j]];
         r[jg[j]] += ggl[j] * x[i];
         j++;
      }
      r[i] = pr[i] - r[i];
   }
   r = Mult(di, N, r);
   z = Mult(di, N, r);
   p = Mult(di, ggl, jg, ig, N, z);
   p = Mult(di, N, p);
   for (i = 0; i != MAXITER && norm_square >= tol; i++)
   {
      r_old = ScalMult(r, r);
      pr_old = ScalMult(p, r);
      p_old = ScalMult(p, p);
      alpha = pr_old / p_old;
      for (j = 0; j < x.size(); j++)
      {
         x[j] = x[j] + alpha * z[j];
         r[j] = r[j] - alpha * p[j];
      }
      temp1 = Mult(di, N, r);
      temp = Mult(di, ggl, jg, ig, N, temp1);
      temp = Mult(di, N, temp);
      beta = -1 * ScalMult(p, temp) / p_old;
      for (j = 0; j < x.size(); j++)
      {
         z[j] = temp1[j] + beta * z[j];
         p[j] = temp[j] + beta * p[j];
      }
      for (j = 0; j < r.size(); j++) {
         norm_square = r_old - pow(alpha, 2) * p_old; // Использовать предыдущую невязку
      }
      cout << "iterations = " << i + 1 << endl;
      cout << "norm_square = " << norm_square << endl;
   }
   if (norm_square >= tol) cout << "exit by max_iterations" << endl;
   else cout << "exit by " << i << " iterrations" << endl;
   /*for (i = 0; i < x.size(); ++i) {
      cout << scientific << setprecision(15) << x[i] << " ";
   }*/
   Write_toFile(x);
}

int main()
{
   setlocale(LC_ALL, "RUS");
   string file_readInf = "readInf.txt";
   int N, MAXITER;
   double tol;
   readInf(file_readInf, N, MAXITER, tol);
   vector<double> pr;
   vector<int> ig;
   vector<int> jg;
   vector<double> ggl;
   vector<double> di;
   Initilize(ig, di, ggl, jg, pr, N);
   //LOS(N, di, ggl, jg, pr, ig, MAXITER, tol);
   LOS_diag(N, di, ggl, jg, pr, ig, MAXITER, tol);

   return 0;
}