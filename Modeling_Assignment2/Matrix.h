#pragma once
class Matrix
{
private:
	double **mat;
	int rows, cols;
public:
	Matrix();
	Matrix(int r, int c);
	int getrows();
	int getcols();
	void add(int i, int j, double v);
	void set(int i, int j, double v);
	double get(int i, int j);
	void printMatrix();
	double** getMatrix();
     void MatrixMultiplication(double **A, int order, double **B, double **result);
	 void MatrixInversion(double **A, int order, double **Y);
	 int GetMinor(double **src, double **dest, int row, int col, int order);
	 double CalcDeterminant(double **mat, int order);
	
	~Matrix();
};

