#pragma once

/*
	This is a Fourier transformer 2D. It allows to work with function on a square [0, 1] * [0, 1]
	with a specific boundary condition:

	f (x=0) = 0
	f'(x=1) = 0
	f'(y=0) = 0
	f (y=1) = 1
*/

#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<vector>
#include<omp.h>


const double PI = 4 * atan(1);


typedef double(*REAL_FUNC)(double, int); //real function pointer
typedef double(*SIMPLE_REAL_FUNC)(double); //real function pointer
typedef double(*REAL_FUNC_2D)(double, double); //2d-real function pointer

//key functions
double func_x(double x, int k);
double func_y(double x, int k);

//scalar product calculator
double scal_prod1D(int n, double h, const double *f, const double *g, bool bias = false);



//to perform furier transform. returns 0 if everything is well
int f2c(int n, double hx, double hy, const double* f, double* c);

//to perform inverced transform. returns 0 if everything is well
int c2f(int n, double hx, double hy, double* f, const double* c);
//this function predicts values of the function in points betwwen the mesh
double* restore_function(int n, double hx, double hy, const double *c, int mult);

//auxilirary functions
void fill_array(int n, double h, double *a, REAL_FUNC f); //fills an array with a function
void fill_array(int n, double h, double* g, SIMPLE_REAL_FUNC f); //fills an array with a function
void fill_array(int n, double hx, double hy, double *a, REAL_FUNC_2D f);

void fill_with_zeros(int n, double* f);
