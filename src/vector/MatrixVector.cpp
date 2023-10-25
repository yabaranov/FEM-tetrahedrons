#include <iostream>
#include "MatrixVector.h"

Vector Vector::operator+(const Vector& v) const
{
	Vector result(v.size());
	for (int i = 0; i < v.size(); i++)
		result[i] = (*this)[i] + v[i];
	return result;
}

Vector Vector::operator-(const Vector& v) const
{
	Vector result(v.size());
	for (int i = 0; i < v.size(); i++)
		result[i] = (*this)[i] - v[i];
	return result;
}

Vector Vector::operator*(double value) const
{
	Vector result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		result[i] = value * (*this)[i];
	return result;
}

Vector Vector::operator/(double value) const
{
	Vector result;
	for (int i = 0; i < (*this).size(); i++)
		result[i] = (*this)[i] / value;
	return result;
}

double Vector::operator*(const Vector& v) const
{
	double result = 0;
	for (int i = 0; i < (*this).size(); i++)
		result += v[i] * (*this)[i];
	return result;
}

Vector Vector::operator-() const
{
	Vector result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		result[i] = -(*this)[i];
	return result;
}

Vector& Vector::operator+=(const Vector& v) 
{
	for (int i = 0; i < (*this).size(); i++)
		(*this)[i] += v[i];
	return *this;
}

Vector& Vector::operator-=(const Vector& v)
{
	for (int i = 0; i < (*this).size(); i++)
		(*this)[i] -= v[i];
	return *this;
}

Vector& Vector::operator=(double value)
{
	for (int i = 0; i < (*this).size(); i++)
		(*this)[i] = value;
	return (*this);
}


Matrix Matrix::operator/(double value) const
{
	Matrix result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			result[i][j] = (*this)[i][j] / value;
	return result;

}

Matrix& Matrix::operator/=(double value)
{
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			(*this)[i][j] /= value;
	return (*this);
}


Matrix& Matrix::operator*=(double value)
{
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			(*this)[i][j] *= value;
	return (*this);
}


Matrix& Matrix::operator=(double value)
{
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			(*this)[i][j] = value;
	return (*this);
}


Vector Matrix::operator*(const Vector& v) const
{
	Vector result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
	{
		result[i] = 0;
		for (int j = 0; j < (*this)[i].size(); j++)
			result[i] += (*this)[i][j] * v[j];
	}
	return result;
}


double Matrix::InverseMatrix()
{
	double eps = 1e-15;
	double detA = 1;
	Matrix inverseA((*this).size());

	auto M = (*this);

	for (int i = 0; i < (*this).size(); i++)
	{
		for (int j = 0; j < i; j++)
			inverseA[i][j] = 0;

		inverseA[i][i] = 1;

		for (int j = i + 1; j < (*this)[i].size(); j++)
			inverseA[i][j] = 0;
	}

	for (int i = 0; i < (*this).size(); i++)
	{
		double max = abs(M[i][i]);
		int maxIndex = i;
		for (int j = i + 1; j < (*this)[i].size(); j++)
		{
			if (abs(M[j][i]) > max)
			{
				max = abs(M[j][i]);
				maxIndex = j;
			}
		}

		if (max < eps)
		{
			std::cout << "Inverse matrix not found! Zero column ";
			std::cout << maxIndex << std::endl;
			return 0;
		}

		std::swap(M[i], M[maxIndex]);
		std::swap(inverseA[i], inverseA[maxIndex]);

		detA *= M[i][i];

		for (int j = 0; j < (*this)[i].size(); j++)
		{
			if (j != i)
			{
				double temp = M[j][i] / M[i][i];
				for (int k = 0; k < (*this)[i].size(); k++)
				{
					M[j][k] -= M[i][k] * temp;
					inverseA[j][k] -= inverseA[i][k] * temp;
				}
			}
		}
	}


	for (int i = 0; i < (*this).size(); i++)
	{
		double temp = M[i][i];
		for (int j = 0; j < (*this)[i].size(); j++)
		{
			inverseA[i][j] /= temp;
		}
	}

	(*this) = inverseA;

	return detA;
}


const Matrix& Matrix::Transpose()
{
	for (int i = 0; i < (*this).size(); i++)
		for (int j = i + 1; j < (*this)[i].size(); j++)
			std::swap((*this)[i][j], (*this)[j][i]);
	return (*this);

}

Matrix Matrix::operator-(const Matrix& M) const
{
	Matrix result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			result[i][j] = (*this)[i][j] - M[i][j];
	return result;
}


Matrix Matrix::operator+(const Matrix& M) const
{
	Matrix result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			result[i][j] = M[i][j] + (*this)[i][j];
	return result;
}


Matrix Matrix::operator*(double value) const
{
	Matrix result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			result[i][j] = value * (*this)[i][j];
	return result;

}

Matrix& Matrix::operator+=(const Matrix& M)
{
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			(*this)[i][j] += M[i][j];
	return (*this);
}

Matrix& Matrix::operator-=(const Matrix& M)
{
	for (int i = 0; i < (*this).size(); i++)
		for (int j = 0; j < (*this)[i].size(); j++)
			(*this)[i][j] -= M[i][j];
	return (*this);
}