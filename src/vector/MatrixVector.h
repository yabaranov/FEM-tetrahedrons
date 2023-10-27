#pragma once

#include <vector>
#include <cmath>

class Vector : public std::vector<double>
{
public:

	Vector(int n = 0) : std::vector<double>(n) {}
	Vector(const std::vector<double>& v) : std::vector<double>(v) {}
	Vector(int n, double value) : std::vector<double>(n, value) {}

	Vector operator+(const Vector& v) const;
	Vector operator-(const Vector& v) const;
	Vector operator*(double value) const;
	Vector operator/(double value) const;
	double operator*(const Vector& v) const;
	Vector operator-() const;
	Vector& operator+=(const Vector& v);
	Vector& operator-=(const Vector& v);
	Vector& operator=(double value);
	friend Vector operator*(double value, const Vector& v) { return v * value; }
	double EuclideanNorm() { return std::sqrt((*this) * (*this)); }
};

class Matrix : public std::vector<Vector>
{
public:
	Matrix(int n = 0) : std::vector<Vector>(n, Vector(n)) {}
	Matrix(const std::vector<Vector>& M) : std::vector<Vector>(M) {}
	Matrix(int n, double value) : std::vector<Vector>(n, Vector(n, value)) {}

	Vector operator*(const Vector& v) const;
	Matrix operator*(double value) const;
	Matrix operator/(double value) const;
	Matrix operator+(const Matrix& M) const;
	Matrix operator-(const Matrix& M) const;
	Matrix& operator=(double value);
	Matrix& operator+=(const Matrix& M);
	Matrix& operator-=(const Matrix& M);
	Matrix& operator*=(double value);
	Matrix& operator/=(double value);
	const Matrix& Transpose();
	double InverseMatrix();
	friend Matrix operator*(double value, const Matrix& M) { return M * value; }
};