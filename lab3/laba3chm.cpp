#include <stdio.h>
#include <tchar.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "slau.h"




int main()
{
	slau sl1;
	//sl1.createm1(10);
	double start,end;
	sl1.readA();
	sl1.readb("pr.txt");

	sl1.createA();
	sl1.printA();
	
	start=clock();
	int k=sl1.calcxlos2();
	end=clock();
	sl1.printx("x.txt");
	printf("\ntime=%lf  ,k=%d",(end-start)/1000.0,k);
	getchar();
	
}