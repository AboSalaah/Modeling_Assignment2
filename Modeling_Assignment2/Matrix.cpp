#include "Matrix.h"
#include <iostream>
using namespace std;

Matrix:: Matrix(){

}
Matrix::Matrix(int r, int c)
{
	rows = r;
	cols = c;

	mat= new double *[r];
	for (int i = 0; i < r; i++)
		mat[i] = new double[c];

	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++) {
			mat[i][j] = 0;
		}
	}

}
int Matrix::getrows() {
	return rows;
}
int Matrix::getcols() {
	return cols;
}
void Matrix::add(int i, int j, double val)
{
	mat[i][j] += val;
}
void Matrix::set(int i, int j, double val)
{
	mat[i][j] = val;
}
double Matrix::get(int i, int j)
{
	return mat[i][j];
}
void Matrix::printMatrix()
{
	for (int i = 1; i < rows; i++)
	{
		for (int j = 1; j < cols; j++)
		{
			cout << this->get(i, j) << " ";
		}
		cout << endl;
	}
	cout << "******************************************************" << endl;
}

double** Matrix::getMatrix() {
	return mat;
}
// matrix inversion
// the result is put in Y
void Matrix::MatrixInversion(double **A, int order, double **Y)
{
	
	// get the determinant of a
	double det = 1.0 / CalcDeterminant(A, order);

	// memory allocation
	double *temp = new double[(order - 1)*(order - 1)];
	double **minor = new double*[order - 1];
	for (int i = 0; i < order - 1; i++)
		minor[i] = temp + (i*(order - 1));

	for (int j = 0; j < order; j++)
	{
		for (int i = 0; i < order; i++)
		{
			// get the co-factor (matrix) of A(j,i)
			GetMinor(A, minor, j, i, order);
			Y[i][j] = det * CalcDeterminant(minor, order - 1);
			if ((i + j) % 2 == 1)
				Y[i][j] = -Y[i][j];
		}
	}
	
	// release memory
	//delete [] minor[0];
	delete[] temp;
	delete[] minor;
}

// calculate the cofactor of element (row,col)
int Matrix::GetMinor(double **src, double **dest, int row, int col, int order)
{
	// indicate which col and row is being copied to dest
	int colCount = 0, rowCount = 0;

	for (int i = 0; i < order; i++)
	{
		if (i != row)
		{
			colCount = 0;
			for (int j = 0; j < order; j++)
			{
				// when j is not the element
				if (j != col)
				{
					dest[rowCount][colCount] = src[i][j];
					colCount++;
				}
			}
			rowCount++;
		}
	}

	return 1;
}

// Calculate the determinant recursively.
double Matrix::CalcDeterminant(double **mat, int order)
{
	// order must be >= 0
	// stop the recursion when matrix is a single element
	if (order == 1)
		return mat[0][0];

	// the determinant value
	double det = 0;

	// allocate the cofactor matrix
	double **minor;
	minor = new double*[order - 1];
	for (int i = 0; i < order - 1; i++)
		minor[i] = new double[order - 1];

	for (int i = 0; i < order; i++)
	{
		// get minor of element (0,i)
		GetMinor(mat, minor, 0, i, order);
		// the recusion is here!

		det += (i % 2 == 1 ? -1.0 : 1.0) * mat[0][i] * CalcDeterminant(minor, order - 1);
		//det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
	}

	// release memory
	for (int i = 0; i < order - 1; i++)
		delete[] minor[i];
	delete[] minor;

	return det;
}


 void Matrix::MatrixMultiplication(double **A, int order, double **B, double **result)
{
	 for (int i = 1; i < order; ++i) {
		 cout << B[i][1] << endl;;
		 
	}
	for (int i = 0; i < order; i++)
	{
		result[i][0] = 0;
		for (int k = 0; k < order; k++)
		{
			result[i][0] = result[i][0] + (A[i][k] * B[k][0]);
		}

	}
}


Matrix::~Matrix()
{
}
