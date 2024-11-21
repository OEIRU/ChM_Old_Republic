#include <stdio.h>
#include <tchar.h>
#include <math.h>
#include <stdlib.h>

class slau
{
private:
	int n;
	double e;
	int maxiter;
	int *ig;
	int *jg;
	double *gg;
	double *di;
	double *b;
	double *x;
	double *L;

	double **A;

	bool isa;
	bool isb;
	bool isx;
public:
	slau();
	int readA();
	int readb(char *filename);

	double multvectors(double *x,double *y);
	void multAvect(double *vect,double *temp);
	double calcnormvect(double *temp);

	void getLU();
	void solveslauLx(double *b,double *x);
	
	int printx(char *filename);

	int calcxlos1();
	int calcxlos2();
	int calcxmsg();

	void createTest(int n);
	void createm1(int n);
	void createm1a(int n);
	void createm2(int n);

	int createA();
	void printA();
};