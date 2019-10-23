#include "fourier.h"
#include<algorithm>
#include<cstdlib>
#include<string>

using namespace std;

double func1(double x, double y) {
	return 2*sin(0.5 * PI * x) * cos(0.5 * PI * y);
}

double func2(double x, double y) {
        return exp(x+y)*(-x)*(y-1);
}

double unit(double x, double y) {
        return 1+x-x+y-y;
}

void test();
double deep_test(int n);

double max_deriv(int n, const double *f, const double* g);
void write_table(int n, const double *a, const string & q);
void take_string(int n, int k, const double *a, double *ans);

int main() {

	int n = 10;


	//check the matrix of scalar products

	double* g1, *g2;
	double *t1, *t2;

	g1 = new double[n*(n+1)];
	g2 = new double[n*(n+1)];

	t1 = new double[n + 1];
	t2 = new double[n + 1];


        double hx = 2. / (2 * n - 1), hy = 1. / n;
	fill_array(n, hx, g1, func_x);
	fill_array(n, hy, g2, func_y);

	//axis x
	ofstream fout("matrix_x.txt");
	

	for (int i = 0; i < n-1; i++) {
		for (int j = 0; j < n-1; j++) {
			
			take_string(n+1, i, g1, t1);
			take_string(n+1, j, g1, t2);
			fout.width(14);
			fout << scal_prod1D(n, hx, t1, t2)<< " ";
		}
		fout << "\n";
	}
	fout.close();
	//axis y
	fout = ofstream("matrix_y.txt");

	for (int i = 0; i < n-1; i++) {
		for (int j = 0; j < n-1; j++) {
			take_string(n+1, i, g2, t1);
			take_string(n+1, j, g2, t2);
			fout.width(14);
			fout << scal_prod1D(n, hy, t1, t2, true)<<" ";
		}
		fout << "\n";
	}
	fout.close();

	delete[] g1; delete[] g2;
        delete[]  t1; delete[] t2;


	

	test();

        ofstream s("deep_research.txt");
        s<<"Research of convergence for the function\n";
        s<<"n            div";
        for(int n = 10; n < 300; n*=2) {
               s<<n;
               s.width(12);
               s<<deep_test(n);
               s<<"\n";
        }

        s.close();

	return 0;

}

void test() {
        int N = 500;
	double hx, hy;
	ofstream fout("research.txt");
	fout << "The research of quality for the function  (0.5 * (x - 1) ^ 2 - 0.5) * (0.25 * y ^ 4 - 0.25) \n";
	fout << "n        err  \n";
	double * f, *c, *tr;

	f = new double[2*N*N];
	c = new double[2*N*N];
	tr = new double[2*N*N];

        fill_with_zeros(2*N*N, c);

	//fill_with_zeros(2 * N*N, c);
	for (int n = 2; n < N; n *= 2) {

		cout << "Testing. Number of segments = " << n << "\n";
		hx = 2. / (2 * n - 1);
		hy = 1. / ( n);

		fill_array(n+1, hx, hy, f, func2);
		cout << "	fourier transform\n";
		f2c(n, hx, hy, f, c);
		cout << " inversed transform\n";
		c2f(n, hx, hy, tr, c);
		fout.width(7); fout << n << " ";
		fout << max_deriv(n+1, f, tr)  << "\n";

		if (n < 4) {
			write_table(n+1, f, "function.txt");
			write_table(n+1, tr, "transform.txt");
			write_table(n, c, "coefs.txt");
		}
	}
	fout.close();
	delete[] f; delete[] c; delete[] tr;
}



double max_deriv(int n, const double* f, const double* g) {
	double c = 0;
	double t;
	for (int i = 0; i < n*n; i++) {
		t = fabs(f[i] - g[i]);
		if (t > c) c = t;
	}
	return c;
}

void write_table(int n, const double *a, const string & q) {
	ofstream f(q);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			f.width(12);
			f << a[n*i + j] << " ";
		}
		f << "\n";
	}
	f.close();
}

void take_string(int n, int k, const double *a, double *ans) {
	for (int i = 0; i < n; i++) ans[i] = a[n*k + i];
}

double deep_test(int n) {
    const int mult = 2;
    double hx = 2./(2*n-1), hy=1./n;
    ofstream fout("research.txt");
    fout << "The research of quality for the function  (0.5 * (x - 1) ^ 2 - 0.5) * (0.25 * y ^ 4 - 0.25) \n";
    fout << "n        err  \n";
    double * f, *c, *tr, *f_ext;

    f =  new double[2*n*n];
    c =  new double[2*n*n];
    f_ext = new double[2*n*n*mult*mult];


    fill_with_zeros(2*n*n, c);
    fill_array(n, hx, hy, f, func2);
    fill_array(n*mult, hx/mult, hy/mult, f_ext, func2);

    f2c(n, hx, hy, f, c);
    tr = restore_function(n, hx, hy, c, mult);

    return max_deriv(n*mult, f_ext, tr);


}
