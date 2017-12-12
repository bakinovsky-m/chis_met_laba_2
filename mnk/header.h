#ifndef HEADER_H
#define HEADER_H

#include <vector>
#include <string>

class Mat{
public:
	Mat() = default;
	Mat(int d);
	Mat(std::vector<std::vector<double>> mat);

	void resize(int d);

	std::vector<double>& operator[](int row);
	std::vector<double> operator[](int row) const;
	Mat operator*(Mat& rhs) const;

	void fillFromFile(std::string str);

	// fiedoubles
	int dim = 0;
	std::vector<double> svobChleny;

	std::vector<double> ownNumbers;
	std::vector<std::vector<double>> ownVectors;

// private:
	std::vector<std::vector<double>> mat;
};

Mat methodGaussJordan(const Mat m);
Mat multRow(const Mat& m, const int ind, const double value);
Mat divRow(const Mat& m, const int ind, const double value);
Mat swapRows(const Mat& m, const int ind1, const int ind2);

void print(const Mat m);
void print(const std::vector<double> vec);

double getDeterminant(const Mat odouble_matr);
Mat getMWithoutIJ(const Mat mas, const int i, const int j);
Mat getMWithIColChangedToSvobChlen(const Mat matr, const int index_stolba);
std::vector<double> kramer(const Mat matr);
std::vector<double> relaxation(const Mat A, const double parameter);
void LUdecompos(const Mat A, Mat &L, Mat &U);

std::vector<std::pair<double, std::vector<double>>> powMethod(const Mat& A, const std::vector<double> svobChleny);
std::vector<double> unmMatNaVec(const Mat& A, const std::vector<double>& vec);
std::vector<double> methodGauss(const Mat& A0, const double l);

Mat umnMat(const Mat & m1, const Mat & m2);

std::vector<double> normal(const std::vector<double> y);
#endif
