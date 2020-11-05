#pragma once
#include<stdio.h>
#include<conio.h>
#include <stdlib.h>
#include <string.h>


int multiplyM(double** A, int NofRowsA, int NofColumnA, double** B, int NofRowsB, int NofColumB, double** outputM);
double determinant(double** A, int order);
void cofactor(double** MIn, double** facOut, int order);
void transpose(double** m, double** m_T, int num_col, int num_raw);
int inverse(double** A, int order, double** InvA);
int polynomialCurvFit( double* x ,  double*  y, int length, int power_of_curve, double* coeff );

//Solve linear system
//double xout is the closet solution of Ax = b
int SolveLinearSys(double** A, double** b, int Num_x_row, int Num_x_col, double* xout);

//this one is slower. Not accurate, but can solve all case. use it when svd method fail
void linearinv(int N, double y[], double **A1, double epsi, double x[], int* ni);

//SVD API
void svdcmp(double **a, int m, int n, double wout[], double **vout);
