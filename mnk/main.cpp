#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "header.h"

using namespace std;

double f(const double x){
	// return pow(x, 3) + x + 5;
	return log(pow(x, 3)) + 3 * pow(x, 2) + 5;
}

double pol(const double x, const double m, const vector<double> c){
	double res = 0;
	for(int i = 0; i <= m; ++i){
		res += c[i] * pow(x, i);
	}
	return res;
}

int main(){
	int n = 100;
	double left = 1;
	double right = 11;
	double step = (right - left)/n;
	double poly_m = 4;

	vector<double> x(n + 1);

	x[0] = left;
	for(int i = 1; i <= n; ++i){
		x[i] = x[i-1] + i*step;
	}

	Mat mat(poly_m+1);
	vector<double> svb(poly_m+1,0);

    for (int i(0); i <= poly_m; ++i) {
 	    for (int j(0); j <= poly_m; ++j) {
	        for (int k(0); k <= n; ++k) {
	    	    mat[i][j] += pow(x[k], i + j);
	        }
	        mat[i][j] /= (n + 1);
	    }
	    for (int k(0); k <= n; ++k) {
	        svb[i] += pow(x[k], i)*f(x[k]);
	    }
	    svb[i] /= (n + 1);
	}

	vector<double> t(0);
	for(int i = 1; i <= mat.dim; ++i){
		t.push_back(i);
	}
	mat.svobChleny = t;
	vector<double> c = relaxation(mat, 1);
	// vector<double> c = methodGauss(mat, 0.5);

	ofstream fil("gnu.plot", ofstream::out);

	for(int i = 0; i <= n; ++i){
		cout << x[i] << " " << f(x[i]) << " " << pol(x[i], poly_m, c) << endl;
		fil << x[i] << " " << f(x[i]) << " " << pol(x[i], poly_m, c) << endl;
	}
	fil.close();

	return 0;
}