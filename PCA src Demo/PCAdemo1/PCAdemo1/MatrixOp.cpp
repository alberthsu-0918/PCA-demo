#include "StdAfx.h"
#include <math.h>
#include "MatrixOp.h"


#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))

//----------Function group: basic matrix operation: mutilpy, add and substract

//  Function multiplyM()
//  Parameters:
//           double** A          : a m*n input matrix
//           int NofRowsA        : the number of rows of matrix A. For a m*n martix, NofRowsA = m 
//           int NofColumnA      : the number of columns of matrix A. For a m*n martix, NofColumA = n 

int multiplyM(double** A, int NofRowsA, int NofColumnA, double** B, int NofRowsB, int NofColumB, double** outputM)
{
	
	int err = -1;
	double sum = 0;
	int row1, column1, row2, column2;

	row1 = NofRowsA;
	column1 = NofColumnA;
	row2 = NofRowsB;
	column2 = NofColumB;

	if(column1 !=  row2)
	{
		printf("A's number of columns not equals to B's number of rows\n");
		err = -1;
	}
	else
	{
		for (int i = 0; i < row1; i++) 
		{
			for (int j = 0; j < column2; j++) 
			{
				for (int k = 0; k < row2; k++) 
				{
					sum = sum + A[i][k]*B[k][j];
				}

				outputM[i][j] = sum;
				sum = 0;
			}
		}

		err = 0;
	}

	return err;
}


//----------Function group: function for inverse a matrix. Includes inverse, transpose, and determine operation
double determinant(double** A, int order)
{
	double s=1;
	double** b;
	int i,j,m,n,c;
	double det;
	
	/*
	double **AA;
	AA = new double* [order];
	for(int i=0;i<order;i++)
	{
		AA[i] = new double[order];
	}
	*/
	if (order==1)
	{
		return (A[0][0]);
	}
	else
	{
		det=0;
		b = new double*[order];
		for(i=0;i<order;i++)
			b[i] = new double[order];

		for (c=0;c<order;c++)
		{
			m=0;
			n=0;
			for (i=0;i<order;i++)
			{
				for (j=0;j<order;j++)
				{
					b[i][j]=0;
					if (i != 0 && j != c)
					{
						b[m][n]=A[i][j];
						if (n<(order-2))
							n++;
						else
						{
							n=0;
							m++;
						}
					}
				}
			}
			det=det + s * (A[0][c] * determinant(b,order-1));
			s=-1 * s;
		}
	}

	for(int i=0;i<order;i++)
		delete [] b[i];
	delete [] b;
	return (det);

}

void transpose(double** m, double** m_T, int num_col, int num_raw)
{
	int i,j, n, c;
	//double b[25][25],inverse[25][25],d;

	n = num_col;
	c = num_raw;
	for (i=0;i<n;i++)
	{
		for (j=0;j<c;j++)
		{
			m_T[i][j]= m[j][i];
		}
	}
}

void cofactor(double** MIn, double** facOut, int order)
{
	double powOrder = -1.000000000;
	double** b;
	b = new double* [order];
	for(int i=0; i<order;i++)
		b[i] = new double[order];

	int p,q,m,n,i,j;
	for (q=0;q<order;q++)
	{
		for (p=0;p<order;p++)
		{
			m=0;
			n=0;
			for (i=0;i<order;i++)
			{
				for (j=0;j<order;j++)
				{
					if (i != q && j != p)
					{
						b[m][n]=MIn[i][j];
						if (n<(order-2))
							n++;
						else
						{
							n=0;
							m++;
						}
					}
				}
			}
			facOut[q][p]=pow(powOrder,q + p) * determinant(b,order-1);
		}
	}

	for(int i=0;i<order;i++)
		delete [] b[i];
	delete [] b;
}

int inverse(double** A, int order, double** InvA)
{
	int err = -1;
	double d;
	int k = order;
	double** fac;

	fac = new double*[order];
	for(int i=0;i<order;i++)
		fac[i] = new double [order];


	d=determinant(A, k) ;
	//d = 5.278355384573969e-035;
	printf("Determinant of the Matrix = %f",d);
	if (d==0)
		printf("\nInverse of Entered Matrix is not possible\n");
	else
		cofactor(A, fac, k);

	transpose(fac,InvA, k, k);
	
	////
	//d=determinant(A,k);
	for (int i=0;i<k;i++)
	{
		for (int j=0;j<k;j++)
		{
			InvA[i][j]=InvA[i][j] / d;
		}
	}
	printf("\n\n\nThe inverse of matrix is : \n");
	
	for(int i=0;i<order;i++)
		delete [] fac[i];

	delete [] fac;

	err = 0;
	return err;
}


//---------------SVD Group---------------
/*******************************************************************************
Singular value decomposition program, svdcmp, from "Numerical Recipes in C"
(Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
and B.P. Flannery
*******************************************************************************/
double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	return v-nl+NR_END;
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

/******************************************************************************/
void svdcmp(double **a, int m, int n, double wout[], double **vout)
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].

Albert Hsu's revision: 2017/05/10
*******************************************************************************/
{
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	// create / 0 0............ 0  \
	//        | 0 a00 a01 ...  a0n |
	//        | 0 a10 a11 ...  a1n |
	//   U =  |..                  |  
	//        |..              ..  |
	//        \ 0 am0 am1 ...  amn / , same as w and v. Pad 0 in the first row and column

	double ** U = new double* [m +1];
	for(int i=0;i<m+1;i++)
	{
		U[i] = new double[n+1];
		memset(U[i] , 0 , sizeof(double)*(n+1));
		if(i>0)
			memcpy(U[i]+1, a[i-1], sizeof(double)*n);
	}

	double *w = new double [n+1];
	memset(w, 0, sizeof(double)*(n+1));

	double **v = new double* [n+1];
	for(int i=0;i<n+1;i++)
	{
		v[i] = new double [n+1];
		memset(v[i], 0, sizeof(double)*(n+1));
	}

	rv1=dvector(1,n);
	g=scale=anorm=0.0; /* Householder reduction to bidiagonal form */
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(U[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					U[k][i] /= scale;
					s += U[k][i]*U[k][i];
				}
				f=U[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				U[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += U[k][i]*U[k][j];
					f=s/h;
					for (k=i;k<=m;k++) U[k][j] += f*U[k][i];
				}
				for (k=i;k<=m;k++) U[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(U[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					U[i][k] /= scale;
					s += U[i][k]*U[i][k];
				}
				f=U[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				U[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=U[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += U[j][k]*U[i][k];
					for (k=l;k<=n;k++) U[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) U[i][k] *= scale;
			}
		}
		anorm = DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) { /* Accumulation of right-hand transformations. */
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) /* Double division to avoid possible underflow. */
					v[j][i]=(U[i][j]/U[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += U[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) { /* Accumulation of left-hand transformations. */
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) U[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += U[k][i]*U[k][j];
				f=(s/U[i][i])*g;
				for (k=i;k<=m;k++) U[k][j] += f*U[k][i];
			}
			for (j=i;j<=m;j++) U[j][i] *= g;
		} else for (j=i;j<=m;j++) U[j][i]=0.0;
		++U[i][i];
	}
	for (k=n;k>=1;k--) { /* Diagonalization of the bidiagonal form. */
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) { /* Test for splitting. */
				nm=l-1; /* Note that rv1[1] is always zero. */
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; /* Cancellation of rv1[l], if l > 1. */
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=U[j][nm];
						z=U[j][i];
						U[j][nm]=y*c+z*s;
						U[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) { /* Convergence. */
				if (z < 0.0) { /* Singular value is made nonnegative. */
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) printf("no convergence in 30 svdcmp iterations");
			x=w[l]; /* Shift from bottom 2-by-2 minor. */
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0; /* Next QR transformation: */
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=U[jj][j];
					z=U[jj][i];
					U[jj][j]=y*c+z*s;
					U[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);


	memcpy(wout, w+1, sizeof(double)*n);
	delete [] w;

	for(int i=0;i<n;i++)
	{
		memcpy(vout[i], v[i+1]+1, sizeof(double)*n);
		delete [] v[i];
	}
	delete [] v;

	for(int i=0;i<m;i++)
	{
		memcpy(a[i], U[i+1]+1, sizeof(double)*n);
		delete [] U[i];
	}
	delete [] U;

}
//---------------End of SVD---------------
//

//   Input parameters are:
//    N:  number of vector elements, and dimension of the coefficient array
//    y:  vector of data points (with dimension N)
//    A1: pointer to first element of the coefficient array
//    epsi: desired relative accuracy (should be chosen not smaller than 
//                                     the relative noise in the data y)  
//
//   Output parameters are:
//    x: solution vector (with dimension N)
//    ni: number of iterations
//
//
//   The function is called in the form    inv(N,y,&A[0][0],epsi,x,ni);
void linearinv(int N, double y[], double **A1, double epsi, double x[], int* ni)  
{

	int i,k;
	//double A[N][N],a[N],xx[N],yy[N],convgc,dev,devmax;
	double** A;
	double *a, *xx, *yy;
	double convgc,dev,devmax;

	a = new double [N];
	xx = new double [N];
	yy = new double [N];
	A = new double* [N];
	for(int i=0;i<N;i++)
	{
		A[i] = new double [N];
	}

	for (i=0; i<N-1; i=i+1)  {
		for (k=0; k<N-1; k=k+1)  {
			A[i][k]=A1[i][k];
		}
	}

	for (i=0; i<N-1; i=i+1)  {
		a[i]=0;
		for (k=0; k<N-1; k=k+1)  {
			if (k!=i) {a[i]=a[i]+fabs(A[i][k]);}
		}
		if (A[i][i]<0) {a[i]=1.1*a[i]+fabs(A[i][i]);}
		A[i][i]=A[i][i]+a[i];
		x[i]=y[i]/A[i][i];
	}

	convgc=1; *ni=0;

	while (convgc >= epsi) {

		for (i=0; i<N-1; i=i+1)  {
			xx[i]=x[i];
		}

		for (i=0; i<N-1; i=i+1)  {
			x[i]=0;
			for (k=0; k<N-1; k=k+1)  {
				if (k!=i) { x[i]=x[i] +xx[k]*A[i][k]; }
			}
			x[i]=(y[i]+a[i]*xx[i]-x[i])/A[i][i] ;
		}


		for (i=0; i<N-1; i=i+1)  {
			yy[i]=0;
			for (k=0; k<N-1; k=k+1)  {
				yy[i]=yy[i]+ x[k]*A[i][k];
			}
			yy[i]=yy[i]- x[i]*a[i];
		}



		devmax=0;
		for (i=1; i<N-1; i=i+1)  {
			if (xx[i]!=0) { dev=fabs(x[i]/xx[i]-1.0); if (dev>devmax) devmax=dev; }
			if (y[i]!=0)  { dev=fabs(yy[i]/y[i]-1.0); if (dev>devmax) devmax=dev; }
		}

		convgc=devmax;
		*ni=*ni+1;

	}

	delete [] a;
	delete [] xx;
	delete [] yy;
	for(int i=0;i<N;i++)
	{
		delete [] A[i];
	}
	delete [] A;

}

//xout is the closet solution of alinear system Ax = b
// Parameters:
//[IN]  A: A of an equation Ax = b
//[IN]  b: b of an equation Ax = b
//[IN]  Num_x_row: number of rows of A. If A is a m*n matrix, Num_x_row = m
//[IN]  Num_x_col: number of colums of A. If A is a m*n matrix, Num_x_col = n
//[Out] xout: estimation of x of Ax = b
int SolveLinearSys(double** A, double** b, int Num_x_row, int Num_x_col, double* xout)
{
	int err = -1;

	//Allocate SVD matrix memory
	double ** U = new double* [Num_x_row];
	for(int i=0;i<Num_x_row;i++)
	{
		U[i] = new double[Num_x_col];
		memcpy(U[i], A[i], sizeof(double)*(Num_x_col));
	}

	double *SW = new double [Num_x_row];
	memset(SW, 0, sizeof(double)*Num_x_row);

	double **V = new double* [Num_x_row];
	for(int i=0;i<Num_x_row;i++)
	{
		V[i] = new double [Num_x_row];
		memset(V[i], 0, sizeof(double)*Num_x_row);
	}
	
	//SVD - singular value decomposition
	svdcmp(U, Num_x_row, Num_x_col, SW, V);
	

	//---start to caculate coeff---
	double ** U_tT = new double* [Num_x_col];
	for(int i=0;i<Num_x_col;i++)
	{
		U_tT[i] = new double[Num_x_row];
	}

	//tempO = U'*y
	transpose(U, U_tT, Num_x_col, Num_x_row);
	double **tempO = new double*[Num_x_col];
	for(int i=0;i<Num_x_col;i++)
	{
		tempO[i] = new double [1];
	}
	err = multiplyM(U_tT, Num_x_col, Num_x_row, b, Num_x_row, 1, tempO);
	for(int i=0;i<Num_x_col;i++)
	{
		tempO[i][0] = tempO[i][0]/SW[i];
	}

	//coeff_SVD = V_t*tempO./diag(SW_t)
	double **coeff_SVD = new double*[Num_x_col];
	for(int i=0;i<Num_x_col;i++)
	{
		coeff_SVD[i] = new double [1];
	}
	err = multiplyM(V, Num_x_col, Num_x_col, tempO, Num_x_col, 1, coeff_SVD);
	//---caculation done---

	//return x
	for(int i=0;i<Num_x_col;i++)
	{
		xout[i] = coeff_SVD[i][0];
	}

	//Clean memory
	for(int i=0;i<Num_x_row;i++)
	{
		delete [] U[i];
		delete [] V[i];
	}
	delete [] U;
	delete [] V;
	delete [] SW;

	for(int i=0;i<Num_x_col;i++)
	{
		delete [] U_tT[i];
		delete [] tempO[i];
		delete [] coeff_SVD[i];
	}
	delete [] U_tT;
	delete [] tempO;
	delete [] coeff_SVD;

	return err;
}

//Function polynomialCurvFit is used to estimate polynomial curve's coefficient, a0~an: Y_est(i) = a0 + a1*x(i) + a2*x(i)^2 + ... + an*x(i)^n
//Parameters:
//(in) x: data array in x axis 
//(in) y: data array in y axis
//(in) length: data length
//(in) power_of_curve: the power of the estimated polynomial curve (=n)
//(out) coeff: the estimation coefficient. coeff =  [a0 a1 a2 ... an]
int  polynomialCurvFit( double* x ,  double*  y, int length, int power_of_curve, double* coeff )
{
	int err = -1;
	double **X,  **Y;
	
	//coeff : from low degree to high degree
    //coeff = inv((X'*X))*X'*y; PS:This way is abondended due to ill-condition of X'*X, which is a close singular matrix(hard to inverse)

	//Allocate memory
	X = new double* [length];	
	Y = new double* [length];

	for(int i=0;i<length; i++)
	{
		X[i] = new double [power_of_curve+1];
		Y[i] = new double [1];
		Y[i][0] = y[i];
	}

	//build Matrix X: where X is a n*m matrix, [x(n)^m], n = 0...length-1; m = 0.....power of curve.
	for(int i=0; i<length; i++)
	{
		X[i][0] = 1;
	}

	for(int i=1; i<power_of_curve+1; i++)
	{
		for(int j=0;j<length;j++)
		{
			X[j][i] = x[j]*X[j][i-1];
		}
	}


	double* xout = new double [power_of_curve+1];
	err =  SolveLinearSys(X, Y, length, power_of_curve+1, coeff);


	// clean memory
	for(int i=0;i<length;i++)
	{
		delete [] X[i];
		delete [] Y[i];
	}

	delete [] X;
	delete [] Y;

	return err;
}

int svd(int m,int n,int withu,int withv,double eps,double tol,
        double *a, double *q, double *u, double *v, double *vt)
{
    
    int i,j,k,l,l1,iter,retval;
    double c,f,g,h,s,x,y,z;
    double *e;
    
    e = (double *)calloc(n,sizeof(double));
    
    retval = 0;

    /* Copy 'a' to 'u' */    
    for (i=0;i<m;i++) {
        for (j=0;j<n;j++)
            u[i*n + j] = a[i*n + j];
    }
    
    /* Householder's reduction to bidiagonal form. */
    g = x = 0.0;
        
    for (i=0;i<n;i++) {
        
        e[i] = g;
        s = 0.0;
        l = i+1;
        
        for (j=i;j<m;j++)
            s += (u[j*n+i]*u[j*n+i]);
        
        if (s < tol)
            g = 0.0;
        
        else {
            
            f = u[i*n+i];
            
            g = (f < 0) ? sqrt(s) : -sqrt(s);
            
            h = f * g - s;
            
            u[i*n+i] = f - g;
            
            for (j=l;j<n;j++) {
                
                s = 0.0;
                
                for (k=i;k<m;k++)
                    s += (u[k*n+i] * u[k*n+j]);
                
                f = s / h;
                
                for (k=i;k<m;k++)
                    u[k*n+j] += (f * u[k*n+i]);
                
            }
            /* end j */
        }
        /* end s */
        q[i] = g;
        
        s = 0.0;
        
        for (j=l;j<n;j++)
            s += (u[i*n+j] * u[i*n+j]);
        
        if (s < tol)
            g = 0.0;
        
        else {
            
            f = u[i*n+i+1];
            
            g = (f < 0) ? sqrt(s) : -sqrt(s);
            
            h = f * g - s;
            
            u[i*n+i+1] = f - g;
            
            for (j=l;j<n;j++) 
                e[j] = u[i*n+j]/h;
            
            for (j=l;j<m;j++) {
                
                s = 0.0;
                
                for (k=l;k<n;k++) 
                    s += (u[j*n+k] * u[i*n+k]);
                
                for (k=l;k<n;k++)
                    u[j*n+k] += (s * e[k]);
                
            }
            /* end j */
        }
        /* end s */
        y = fabs(q[i]) + fabs(e[i]);
                                 
        if (y > x)
            x = y;
        
    }
    /* end i */

    /* accumulation of right-hand transformations */
    if (withv) {
        
        for (i=n-1;i>=0;i--) {
            
            if (g != 0.0) {
                
                h = u[i*n+i+1] * g;
                
                for (j=l;j<n;j++)
                    v[j*n+i] = u[i*n+j]/h;
                
                for (j=l;j<n;j++) {
                    
                    s = 0.0;
                    
                    for (k=l;k<n;k++) 
                        s += (u[i*n+k] * v[k*n+j]);
                    
                    for (k=l;k<n;k++)
                        v[k*n+j] += (s * v[k*n+i]);
                    
                    
                }
                /* end j */
            }
            /* end g */
            for (j=l;j<n;j++)
                v[i*n+j] = v[j*n+i] = 0.0;
            
            v[i*n+i] = 1.0;
            
            g = e[i];
            
            l = i;
            
        }
        /* end i */
        
    }
    /* end withv, parens added for clarity */

    /* accumulation of left-hand transformations */
    if (withu) {
        
        for (i=n;i<m;i++) {
            
            for (j=n;j<m;j++)
                u[i*n+j] = 0.0;
            
            u[i*n+i] = 1.0;
            
        }
        
    }
    
    if (withu) {
        
        for (i=n-1;i>=0;i--) {
            
            l = i + 1;
            
            g = q[i];
            
            for (j=l;j<m;j++)  /* upper limit was 'n' */
                u[i*n+j] = 0.0;
            
            if (g != 0.0) {
                
                h = u[i*n+i] * g;
                
                for (j=l;j<m;j++) {
                    /* upper limit was 'n' */
                    s = 0.0;
                    
                    for (k=l;k<m;k++)
                        s += (u[k*n+i] * u[k*n+j]);
                    
                    f = s / h;
                    
                    for (k=i;k<m;k++) 
                        u[k*n+j] += (f * u[k*n+i]);
                    
                }
                /* end j */
                for (j=i;j<m;j++) 
                    u[j*n+i] /= g;
                
            }
            /* end g */
            else {
                
                for (j=i;j<m;j++)
                    u[j*n+i] = 0.0;
                
            }
            
            u[i*n+i] += 1.0;
            
        }
        /* end i*/
    }
    /* end withu, parens added for clarity */

    /* diagonalization of the bidiagonal form */
    eps *= x;
    
    for (k=n-1;k>=0;k--) {
        
        iter = 0;
        
    test_f_splitting:
        for (l=k;l>=0;l--) {
            
            if (fabs(e[l]) <= eps) goto test_f_convergence;
            
            if (fabs(q[l-1]) <= eps) goto cancellation;
            
        }
        /* end l */

        /* cancellation of e[l] if l > 0 */
    cancellation:
        c = 0.0;
        
        s = 1.0;
        
        l1 = l - 1;
        
        for (i=l;i<=k;i++) {
            
            f = s * e[i];
            
            e[i] *= c;
            
            if (fabs(f) <= eps) goto test_f_convergence;
            
            g = q[i];
            
            h = q[i] = sqrt(f*f + g*g);
            
            c = g / h;
            
            s = -f / h;
            
            if (withu) {
                
                for (j=0;j<m;j++) {
                    
                    y = u[j*n+l1];
                    
                    z = u[j*n+i];
                    
                    u[j*n+l1] = y * c + z * s;
                    
                    u[j*n+i] = -y * s + z * c;
                    
                }
                /* end j */
            }
            /* end withu, parens added for clarity */
        }
        /* end i */
    test_f_convergence:
        z = q[k];
        
        if (l == k) goto convergence;
        

        /* shift from bottom 2x2 minor */
        iter++;
        
        if (iter > 30) {
            
            retval = k;
            
            break;
            
        }
        
        x = q[l];
        
        y = q[k-1];
        
        g = e[k-1];
        
        h = e[k];
        
        f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y);
        
        g = sqrt(f*f + 1.0);
        
        f = ((x-z)*(x+z) + h*(y/((f<0)?(f-g):(f+g))-h))/x;
        
        /* next QR transformation */
        c = s = 1.0;
        
        for (i=l+1;i<=k;i++) {
            
            g = e[i];
            
            y = q[i];
            
            h = s * g;
            
            g *= c;
            
            e[i-1] = z = sqrt(f*f+h*h);
            
            c = f / z;
            
            s = h / z;
            
            f = x * c + g * s;
            
            g = -x * s + g * c;
            
            h = y * s;
            
            y *= c;
            
            if (withv) {
                
                for (j=0;j<n;j++) {
                    
                    x = v[j*n+i-1];
                    
                    z = v[j*n+i];
                    
                    v[j*n+i-1] = x * c + z * s;
                    
                    v[j*n+i] = -x * s + z * c;
                    
                }
                /* end j */
            }
            /* end withv, parens added for clarity */
            q[i-1] = z = sqrt(f*f + h*h);
            
            c = f/z;
            
            s = h/z;
            
            f = c * g + s * y;
            
            x = -s * g + c * y;
            
            if (withu) {
                
                for (j=0;j<m;j++) {
                    
                    y = u[j*n+i-1];
                    
                    z = u[j*n+i];
                    
                    u[j*n+i-1] = y * c + z * s;
                    
                    u[j*n+i] = -y * s + z * c;
                    
                }
                /* end j */
            }
            /* end withu, parens added for clarity */
        }
        /* end i */
        e[l] = 0.0;
        
        e[k] = f;
        
        q[k] = x;
        
        goto test_f_splitting;
        
    convergence:
        if (z < 0.0) {
            
            /* q[k] is made non-negative */
            q[k] = - z;
            
            if (withv) {
                
                for (j=0;j<n;j++)
                    v[j*n+k] = -v[j*n+k];
                
            }
            /* end withv, parens added for clarity */
        }
        /* end z */
    }
    /* end k */
    
    free(e);

    //matrix_transpose(m, n, v, vt);

    return retval;
    
}