#include "slau.h"

slau::slau()
{
	n=0;
	maxiter=0;
	e=1000;
	ig=NULL;
	jg=NULL;
	gg=NULL;
	di=NULL;
	b=NULL;
	x=NULL;
	isx=false;
	isb=false;
	isa=false;
}

int slau::readA()
{
	FILE *f0=fopen("kuslau.txt","r");
	FILE *f1=fopen("ig.txt","r");
	FILE *f2=fopen("gg.txt","r");
	FILE *f3=fopen("jg.txt","r");
	FILE *f4=fopen("di.txt","r");
	if(f0==NULL||f1==NULL||f2==NULL||f3==NULL||f4==NULL)return 0;
	
	if(fscanf(f0,"%d  %d  %lf",&n,&maxiter,&e)==EOF)return 0;

	
	ig=new int[n+1];
	di=new double[n];

	for(int i=0;i<n+1;i++)
	{
		if(fscanf(f1,"%d",&ig[i])==EOF)return 0;
		ig[i]--;
	}

	gg=new double[ig[n]];
	for(int i=0;i<ig[n];i++)
		if(fscanf(f2,"%lf",&gg[i])==EOF)return 0;
	
	jg=new int[ig[n]];
	for(int i=0;i<ig[n];i++)
	{
		if(fscanf(f3,"%d",&jg[i])==EOF)return 0;
		jg[i]--;
	}

	for(int i=0;i<n;i++)
		if(fscanf(f4,"%lf",&di[i])==EOF)return 0;

	x=new double[n];

	for(int i=0;i<n;i++)
		x[i]=0;

	fclose(f0);
	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);
	isa=true;
	return 1;
}


int slau::readb(char *filename)
{
	FILE *f=fopen(filename,"r");
	if(f==NULL)return 0;

	b=new double[n];

	for(int i=0;i<n;i++)
		if(fscanf(f,"%lf",&b[i])==EOF)return 0;

	isb=true;
	return 1;
}


int slau::printx(char *filename)
{
	if(!isx)return 0;
	FILE *f=fopen(filename,"w");
	for(int i=0;i<n;i++)
		fprintf(f,"\n%lf",x[i]);
	return 1;
}

double slau::multvectors(double *x,double *y)
{
	double s=0;
	for(int i=0;i<n;i++)
	{
		s+=x[i]*y[i];
	}
	return s;
}

double slau::calcnormvect(double *temp)
{
	return sqrt(multvectors(temp,temp));
}

void slau::multAvect(double *vect,double *temp)
{
	if(!isa)return;
	int j;
	for(int i=0;i<n;i++)
	{
		temp[i]=di[i]*vect[i];
		for(int k=ig[i];k<ig[i+1];k++)
		{
			j=jg[k];
			temp[i]+=gg[k]*vect[j];
			temp[j]+=gg[k]*vect[i];
		}		
	}
}


void slau::getLU()
{
	L=new double[n];
	for(int i=0;i<n;i++)
		L[i]=sqrt(di[i]);
}

void slau::solveslauLx(double *b,double *x)
{
	for(int i=0;i<n;i++)
		x[i]=b[i]/L[i];
}

int slau::calcxlos1()
{
	if(!isa||!isb)return 0;

	int k;
	double *r,*p,*z;
	r=new double[n];
	p=new double[n];
	z=new double[n];

	double *temp=new double[n];
	multAvect(x,temp);

	for(int i=0;i<n;i++)
	{
		r[i]=b[i]-temp[i];
		z[i]=r[i];
	}

	multAvect(z,temp);

	for(int i=0;i<n;i++)
		p[i]=temp[i];

	double ak,bk;

	for(k=1;k<=maxiter&&multvectors(r,r)>e;k++)
	{	
		ak=multvectors(p,r)/multvectors(p,p);
		for(int i=0;i<n;i++)
		{
			x[i]+=ak*z[i];
			r[i]+=-ak*p[i];
		}

		multAvect(r,temp);

		bk=-multvectors(p,temp)/multvectors(p,p);

		for(int i=0;i<n;i++)
		{
			z[i]=r[i]+bk*z[i];
			p[i]=temp[i]+bk*p[i];
		}
		printf("\n %lf   %lf    otn_nev=%e",multvectors(r,r),e,calcnormvect(r)/calcnormvect(b));
	}

	isx=true;
	return k-1;
}


int slau::calcxmsg()
{
	if(!isa||!isb)return 0;
	int k;
	double *r,*z,*R;
	r=new double[n];
	z=new double[n];
	R=new double[n];

	double *temp=new double[n];
	multAvect(x,temp);

	for(int i=0;i<n;i++)
	{
		r[i]=b[i]-temp[i];
		z[i]=r[i];
	}
	double ak,bk;
	for(k=1;k<=maxiter&&calcnormvect(r)/calcnormvect(b)>e;k++)//multvectors(r,r)>e;k++)
	{
		multAvect(z,temp);
		ak=multvectors(r,r)/multvectors(temp,z);
		
		for(int i=0;i<n;i++)
		{
			x[i]=x[i]+ak*z[i];
			R[i]=r[i]-ak*temp[i];
		}
		
		bk=multvectors(R,R)/multvectors(r,r);

		for(int i=0;i<n;i++)
		{
			z[i]=R[i]+bk*z[i];
			r[i]=R[i];
		}
	}	
	isx=true;
	return k-1;
}


int slau::calcxlos2()
{
	int k;
	if(!isa||!isb)return 0;
	double *r,*z,*p;
	r=new double[n];
	z=new double[n];
	p=new double[n];
	double ak,bk;

	double *temp=new double[n];
	double *temp2=new double[n];
	
	multAvect(x,temp);
	for(int i=0;i<n;i++)
	{
		temp[i]=-temp[i]+b[i];
	}
	getLU();
	solveslauLx(temp,r); 

	solveslauLx(r,z);	

	multAvect(z,temp);
	solveslauLx(temp,p); 

	for(k=1;k<=maxiter&&multvectors(r,r)>e;k++)
	{
		ak=multvectors(p,r)/multvectors(p,p);

		for(int i=0;i<n;i++)
		{
			x[i]=x[i]+ak*z[i];
			r[i]=r[i]-ak*p[i];
		}

		solveslauLx(r,temp);
		multAvect(temp,temp2);
		solveslauLx(temp2,temp);
		bk=-multvectors(p,temp)/multvectors(p,p);

		solveslauLx(r,temp2);
		for(int i=0;i<n;i++)
		{
			z[i]=temp2[i]+bk*z[i];
			p[i]=temp[i]+bk*p[i];
		}
		printf("\n %e   %e    otn_nev=%e",multvectors(r,r),e,calcnormvect(r)/calcnormvect(b));
	}
	
	isx=true;
	return k-1;
}


void slau::createTest(int n)
{
	this->n=n;
	di=new double[n];
	x=new double[n];
	b=new double[n];

	for(int i=0;i<n;i++)
		di[i]=100;
	
	ig=new int[n+1];
	ig[0]=0;
	for(int i=1;i<n+1;i++)
		ig[i]=ig[i-1]+i-1;

	jg=new int[ig[n]];
	jg[0]=0;
	for(int i=1,k=0;i<n;i++)
	{
		for(int k1=0,ikol=i;ikol>0;ikol--,k++,k1++)
			jg[k]=k1;			
	}

	gg=new double[ig[n]];
	for(int i=0;i<ig[n];i++)
	{
		gg[i]=1;
	}	
	
	maxiter=10;
	e=0.00001;
	isa=true;	
	for(int i=0;i<n;i++)
		x[i]=i+1;
	
	multAvect(x,b);
	for(int i=0;i<n;i++)
		x[i]=0;
	isb=true;
}	


void slau::createm1(int n)
{
	int q,ikol,k,k1,j;
	this->n=n;
	di=new double[n];
	x=new double[n];
	b=new double[n];

	ig=new int[n+1];
	ig[0]=0;
	for(int i=1;i<n+1;i++)
		ig[i]=ig[i-1]+i-1;

	jg=new int[ig[n]];
	jg[0]=0;
	for(int i=1,k=0;i<n;i++)
	{
		for(k1=0,ikol=i;ikol>0;ikol--,k++,k1++)
			jg[k]=k1;			
	}

	gg=new double[ig[n]];
		for(int i=0;i<ig[n];i++)
		{
			q=rand()%100;
			if(q<20)gg[i]=0;
			else
				if(q<40)gg[i]=-1;
				else
					if(q<60)gg[i]=-2;
					else
						if(q<80)gg[i]=-3;
						else
							gg[i]=-4;
		}
		
	di=new double[n];
	for(int i=0;i<n;i++)
		di[i]=0;
	for(int i=0;i<n;i++)
	{
		di[i]=0;
		for(k=ig[i];k<ig[i+1];k++)
		{
			di[i]-=gg[k];
		}
	}
	for(int i=0;i<n;i++)
	{
		for(k=ig[i];k<ig[i+1];k++)
		{
			j=jg[k];
			di[j]-=gg[k];
		}
	}

	di[0]+=1;
	maxiter=1000;
	e=1e-40;
	isa=true;

	for(int i=0;i<n;i++)
		x[i]=i+1;
	
	multAvect(x,b);
	for(int i=0;i<n;i++)
		x[i]=0;
	isb=true;
}


void slau::createm1a(int n)
{
	int q,ikol,k,k1,j;
	this->n=n;
	di=new double[n];
	x=new double[n];
	b=new double[n];

	ig=new int[n+1];
	ig[0]=0;
	for(int i=1;i<n+1;i++)
		ig[i]=ig[i-1]+i-1;

	jg=new int[ig[n]];
	jg[0]=0;
	for(int i=1,k=0;i<n;i++)
	{
		for(k1=0,ikol=i;ikol>0;ikol--,k++,k1++)
			jg[k]=k1;			
	}

	gg=new double[ig[n]];
		for(int i=0;i<ig[n];i++)
		{
			q=rand()%100;
			if(q<20)gg[i]=0;
			else
				if(q<40)gg[i]=1;
				else
					if(q<60)gg[i]=2;
					else
						if(q<80)gg[i]=3;
						else
							gg[i]=4;
		}
		
	di=new double[n];
	for(int i=0;i<n;i++)
	{
		di[i]=0;
		for(k=ig[i];k<ig[i+1];k++)
		{
			di[i]+=gg[k];
		}
	}
	for(int i=0;i<n;i++)
	{
		for(k=ig[i];k<ig[i+1];k++)
		{
			j=jg[k];
			di[j]+=gg[k];
		}
	}

	di[0]+=1;
	maxiter=1000;
	e=1e-7;
	isa=true;

	for(int i=0;i<n;i++)
		x[i]=i+1;
	
	multAvect(x,b);
	for(int i=0;i<n;i++)
		x[i]=0;
	isb=true;
}



void slau::createm2(int n)
{
	int q,ikol,k,k1,l,j;
	this->n=n;
	x=new double[n];
	b=new double[n];

	ig=new int[n+1];
	ig[0]=0;
	for(int i=1;i<n+1;i++)
		ig[i]=ig[i-1]+i-1;

	jg=new int[ig[n]];
	jg[0]=0;
	for(int i=1,k=0;i<n;i++)
	{
		for(k1=0,ikol=i;ikol>0;ikol--,k++,k1++)
			jg[k]=k1;			
	}

	gg=new double[ig[n]];
	for(int i=0,k=0;i<n;i++)
	{
		for(k=ig[i];k<ig[i+1];k++)
		{
			j=jg[k];
			gg[k]=1.0/(i+j+1.0);
		}
	}
	
	di=new double[n];
	for(int i=0;i<n;i++)
		di[i]=1.0/(2*i+1.0);

	e=0.001;
	maxiter=1000;
	isa=true;

	for(int i=0;i<n;i++)
		x[i]=i+1;
	
	multAvect(x,b);
	for(int i=0;i<n;i++)
		x[i]=0;
	isb=true;

}


int slau::createA()
{
	if(!isa)return 0;
	A=new double*[n];
	for(int i=0;i<n;i++)
		A[i]=new double[n];

	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			A[i][j]=0;

	int j;
	for(int i=0,k=0;i<n;i++)
	{
		for(k=ig[i];k<ig[i+1];k++)
		{
			j=jg[k];
			A[i][j]=gg[k];
		}
	}

	for(int i=0;i<n;i++)
		for(j=i+1;j<n;j++)
			A[i][j]=A[j][i];

	for(int i=0;i<n;i++)
		A[i][i]=di[i];
	return 1;
}


void slau::printA()
{
	for(int i=0;i<n;i++)
	{
		printf("\n");
		for(int j=0;j<n;j++)
			printf("%lf  ",A[i][j]);
	}
	printf("\n");
	printf("\n");
	for(int i=0;i<ig[n];i++)
		printf("%lf ",gg[i]);

	printf("\n");
	for(int i=0;i<n;i++)
		printf("%d  ",ig[i]);

	printf("\n");
		for(int i=0;i<ig[n];i++)
		printf("%d  ",jg[i]);
}