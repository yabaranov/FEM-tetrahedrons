#ifndef FINITE_MATRIX_VECTOR
#define FINITE_MATRIX_VECTOR
#include <array>
#include <iostream>
template <typename int n>
class FiniteVector : public std::array<double, n>
{

public:
	FiniteVector(const std::array<double, n>& v) : std::array<double, n>(v) {}
	FiniteVector() {}
	~FiniteVector() {}

	FiniteVector<n>& operator=(double value);
	FiniteVector<n> operator+(const FiniteVector<n>& v) const;
	FiniteVector<n> operator-(const FiniteVector<n>& v) const;
	FiniteVector<n> operator*(double value) const;
	FiniteVector<n> operator/(double value) const;
	FiniteVector<n>& operator/=(double value);
	FiniteVector<n>& operator+=(const FiniteVector<n>& v);
	FiniteVector<n>& operator-=(const FiniteVector<n>& v);
	double operator*(const FiniteVector<n>& v) const;
	FiniteVector<n> operator-() const;
	friend FiniteVector<n> operator*(double value, const FiniteVector<n>& v) { return v * value; }
	double distance(const FiniteVector<n>& v) const;
	double EuclideanNorm();
};

template <typename int n>
FiniteVector<n> FiniteVector<n>::operator/(double value) const
{
	FiniteVector<n> result;
	for (int i = 0; i < n; i++)
		result[i] = (*this)[i] / value;
	return result;
}

template <typename int n>
double FiniteVector<n>::EuclideanNorm()
{
	return std::sqrt((*this) * (*this));
}


template <typename int n>
FiniteVector<n>& FiniteVector<n>::operator-=(const FiniteVector<n>& v)
{
	for (int i = 0; i < n; i++)
		(*this)[i] -= v[i];
	return (*this);
}

template <typename int n>
FiniteVector<n>& FiniteVector<n>::operator+=(const FiniteVector<n>& v)
{
	for (int i = 0; i < n; i++)
			(*this)[i]+= v[i];
	return (*this);
}

template <typename int n>
double FiniteVector<n>::distance(const FiniteVector<n>& v) const
{
	double dist = 0;
	for (int i = 0; i < (*this).size(); i++)
		dist += ((*this)[i] - v[i]) * ((*this)[i] - v[i]);
	return std::sqrt(dist);
}

template <typename int n>
FiniteVector<n>& FiniteVector<n>::operator=(double value)
{
	for (int i = 0; i < n; i++)
		(*this)[i] = value;
	return (*this);
}

template <typename int n>
FiniteVector<n> FiniteVector<n>::operator-() const
{
	FiniteVector<n> result;
	for (int i = 0; i < n; i++)
		result[i] = -(*this)[i];

	return result;
}

template <typename int n>
FiniteVector<n>& FiniteVector<n>::operator/=(double value)
{
	for (int i = 0; i < n; i++)
		(*this)[i] /= value;
	return (*this);
}

template <typename int n>
double FiniteVector<n>::operator*(const FiniteVector<n>& v) const
{
	double result = 0;
	for (int i = 0; i < n; i++)
		result += v[i] * (*this)[i];
	return result;
}

template <typename int n>
FiniteVector<n> FiniteVector<n>::operator+(const FiniteVector<n>& v) const
{
	FiniteVector<n> result;
	for (int i = 0; i < n; i++)
		result[i] = (*this)[i] + v[i];
	return result;
}

template <typename int n>
FiniteVector<n> FiniteVector<n>::operator-(const FiniteVector<n>& v) const
{
	FiniteVector<n> result;
	for (int i = 0; i < n; i++)
		result[i] = (*this)[i] - v[i];
	return result;
}

template <typename int n>
FiniteVector<n> FiniteVector<n>::operator*(double value) const
{
	FiniteVector<n> result;
	for (int i = 0; i < n; i++)
		result[i] = value * (*this)[i];
	return result;
}


template <typename int n>
class FiniteMatrix : public std::array<FiniteVector<n>, n>
{

public:
	FiniteMatrix(const std::array<FiniteVector<n>, n>& M) : std::array<FiniteVector<n>, n>(M) {}
	FiniteMatrix() {}
	~FiniteMatrix() {}

	FiniteVector<n> operator*(const FiniteVector<n>& v) const;
	FiniteMatrix<n> operator*(double value) const;
	FiniteMatrix<n> operator/(double value) const;
	FiniteMatrix<n> operator+(const FiniteMatrix<n>& M) const;
	FiniteMatrix<n> operator-(const FiniteMatrix<n>& M) const;
	FiniteMatrix<n>& operator=(double value);
	FiniteMatrix<n>& operator+=(const FiniteMatrix<n>& M);
	FiniteMatrix<n>& operator-=(const FiniteMatrix<n>& M);
	FiniteMatrix<n>& operator*=(double value);
	FiniteMatrix<n>& operator/=(double value);
	const FiniteMatrix<n>& Transpose();
	double InverseMatrix();
	friend FiniteMatrix<n> operator*(double value, const FiniteMatrix<n>& M) { return M * value; }
};

template <typename int n>
FiniteMatrix<n> FiniteMatrix<n>::operator/(double value) const
{
	FiniteMatrix<n> result;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			result[i][j] = (*this)[i][j] / value;
	return result;

}
template <typename int n>
FiniteMatrix<n>& FiniteMatrix<n>::operator/=(double value)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			(*this)[i][j] /= value;
	return (*this);
}

template <typename int n>
FiniteMatrix<n>& FiniteMatrix<n>::operator*=(double value)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			(*this)[i][j] *= value;
	return (*this);
}

template <typename int n>
FiniteMatrix<n>& FiniteMatrix<n>::operator=(double value)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			(*this)[i][j] = value;
	return (*this);
}

template <typename int n>
FiniteVector<n> FiniteMatrix<n>::operator*(const FiniteVector<n>& v) const
{
	FiniteVector<n> result;
	for (int i = 0; i < n; i++)
	{
		result[i] = 0;
		for (int j = 0; j < n; j++)
			result[i] += (*this)[i][j] * v[j];
	}
	return result;
}

template <typename int n>
double FiniteMatrix<n>::InverseMatrix()
{
	double eps = 1e-15;
	double detA = 1;
	FiniteMatrix<n> inverseA;

	auto M = (*this);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
			inverseA[i][j] = 0;

		inverseA[i][i] = 1;

		for (int j = i + 1; j < n; j++)
			inverseA[i][j] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		double max = abs(M[i][i]);
		int maxIndex = i;
		for (int j = i + 1; j < n; j++)
		{
			if (abs(M[j][i]) > max)
			{
				max = abs(M[j][i]);
				maxIndex = j;
			}
		}

		if (max < eps)
		{
			std::cerr << "Inverse matrix not found! Zero column ";
			std::cerr << maxIndex << std::endl;
			return 0;
		}

		std::swap(M[i], M[maxIndex]);
		std::swap(inverseA[i], inverseA[maxIndex]);

		detA *= M[i][i];

		for (int j = 0; j < n; j++)
		{
			if (j != i)
			{
				double temp = M[j][i] / M[i][i];
				for (int k = 0; k < n; k++)
				{
					M[j][k] -= M[i][k] * temp;
					inverseA[j][k] -= inverseA[i][k] * temp;
				}
			}
		}
	}


	for (int i = 0; i < n; i++)
	{
		double temp = M[i][i];
		for (int j = 0; j < n; j++)
		{
			inverseA[i][j] /= temp;
		}
	}

	(*this) = inverseA;

	return detA;
}

template <typename int n>
const FiniteMatrix<n>& FiniteMatrix<n>::Transpose()
{
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			std::swap((*this)[i][j], (*this)[j][i]);
	return (*this);

}
template <typename int n>
FiniteMatrix<n> FiniteMatrix<n>::operator-(const FiniteMatrix<n>& M) const
{
	FiniteMatrix<n> result;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			result[i][j] = (*this)[i][j] - M[i][j];
	return result;
}

template <typename int n>
FiniteMatrix<n> FiniteMatrix<n>::operator+(const FiniteMatrix<n>& M) const
{
	FiniteMatrix<n> result;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			result[i][j] = M[i][j] + (*this)[i][j];
	return result;
}

template <typename int n>
FiniteMatrix<n> FiniteMatrix<n>::operator*(double value) const
{
	FiniteMatrix<n> result;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			result[i][j] = value * (*this)[i][j];
	return result;

}

template <typename int n>
FiniteMatrix<n>& FiniteMatrix<n>::operator+=(const FiniteMatrix<n>& M)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			(*this)[i][j] += M[i][j];
	return (*this);
}

template <typename int n>
FiniteMatrix<n>& FiniteMatrix<n>::operator-=(const FiniteMatrix<n>& M)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			(*this)[i][j] -= M[i][j];
	return (*this);
}

#endif