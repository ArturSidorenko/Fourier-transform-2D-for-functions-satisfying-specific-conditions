/*
	The source file with a 2D transformer
*/

#include "./fourier.h"

const double sq2 = sqrt(2);

double func_x(double x, int k) {
	return sin(x * PI * (k + 0.5)) * sq2;
}

double func_y(double x, int k) {
	return cos(x * PI * (k + 0.5)) * sq2;
}

double scal_prod1D(int n, double h, const double *f,  const double *g, bool bias) {
	double sum = 0;

	if (!bias) {
		for (int i = 0; i < n; i++) sum += f[i] * g[i] * h;
	}
	else {
		for (int i = 1; i < n; i++) sum += f[i] * g[i] * h;
		sum += (f[0] * g[0] + f[n] * g[n])*h*0.5;
	}
	return sum;
}


int f2c(int n, double hx, double hy, const double *f, double *c) {
	double *gx, *gy, *temp;
	gx = new double[n*(n + 1)];
	gy = new double[n*(n+1)];
	temp = new double[n + 1];

	fill_array(n, hx, gx, func_x);
	fill_array(n, hy, gy, func_y);

	for (int i = 0; i < n; i++) {
		#pragma parallel for
		for (int j = 0; j < n+1; j++) temp[j] = scal_prod1D(n, hy, f + (n + 1)*j, gy + (n+1)*i, true);
		#pragma parallel for
		for (int j = 0; j < n-1; j++) c[n*j + i] = scal_prod1D(n, hx, temp, gx + (n + 1)*j);
	}


	delete[] gx;
	delete[] gy;
        delete[] temp;
	return 0;
}

int c2f(int n, double hx, double hy, double *f, const double *c) {
	double *gx, *gy;
	gx = new double[n*n + n];
	gy = new double[n*n + n];

	fill_array(n, hx, gx, func_x);
	fill_array(n, hy, gy, func_y);
	fill_with_zeros((n + 1)*(n + 1), f);

	//a sophisticated formula
	#pragma parallel for
	for (int i = 0; i < n - 1; i++) {
		double temp = 0;
		for (int l = 0; l < n + 1; l++) {
			temp = 0;
			for (int j = 0; j < n; j++) {
				temp += c[i*(n) + j] * gy[(n + 1)*j + l];
			}

			for (int k = 0; k < n+1; k++) {
				f[(n + 1)*k + l] += temp * gx[(n + 1)*i + k];
			}
		}
	}



	delete[] gx;
	delete[] gy;
	return 0;
}


void fill_array(int n, double h, double *a, REAL_FUNC f) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n + 1; j++) {
			a[(n + 1)*i + j] = f(j*h, i);
		}
	}
}

void fill_with_zeros(int n, double *a) {
	for (int i = 0; i < n; i++) a[i] = 0;
}

void fill_array(int n, double hx, double hy, double *a, REAL_FUNC_2D f) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= n; j++) a[(n+1)*i + j] = f(i*hx, j*hy);
	}
	//continuation
	for (int j = 0; j <= n; j++) a[(n)*(n+1) + j] = a[(n-1)*(n+1) + j];
}

void fill_array(int n, double h, double* g, SIMPLE_REAL_FUNC f)
{
	for (int i = 0; i < n; i++) g[i] = f(i*h);
}


//restores function inside the cells of the mesh
double* restore_function(int n, double hx, double hy, const double *c, int mult) {
    //double temp = 0;
	double hx_new = hx / mult;
	double hy_new = hy / mult;

    double *ans;
    ans = new double[(n*mult+1)*(n*mult+1)];

    fill_with_zeros((n*mult + 1)*(n*mult + 1), ans);

    //a sophisticated formula

	#pragma parallel for
    for (int i = 0; i < n - 1; i++) {
		double temp = 0;
            for (int l = 0; l <= n*mult; l++) {
                    temp = 0;
                    for (int j = 0; j < n; j++) {
                            temp += c[i*(n) + j] * func_y(l*hy_new, j);
                    }

                    for (int k = 0; k <= n*mult; k++) {
                            ans[(n*mult+1)*k + l] += temp * func_x(k*hx_new, i);
                    }
            }
    }

    return ans;
}







