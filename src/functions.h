#ifndef FUNCTIONS
#define FUNCTIONS

#include <vector>
#include <functional>
#include <array>
#include "vector/FiniteMatrixVector.h"

enum DIMENSIONS
{
	SIZE_NODE = 3,
	SIZE_ELEMENT = 4,
	SIZE_SUBDOMAIN = 6,
	SIZE_EDGE = 3,
	NUMBER_NODES_CUBE = 8
};

template <typename T>
inline int binarySearch(const std::vector<T>& values, const T& value, int l, int r)
{
	while (l != r)
	{
		int mid = (l + r) / 2 + 1;
		(values[mid] > value) ? r = mid - 1: l = mid;
	}

	return (values[l] == value) ? l : -1;
}

template <typename T, int n>
inline int binarySearch(const std::array<T, n>& values, const T& value, int l, int r)
{
	while (l != r)
	{
		int mid = (l + r) / 2 + 1;
		(values[mid] > value) ? r = mid - 1 : l = mid;
	}

	return (values[l] == value) ? l : -1;
}

inline double derivativeX(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1)
{
	std::array<double, SIZE_NODE> xPlush = args;
	xPlush[0] += h;
	std::array<double, SIZE_NODE> xMinush = args;
	xMinush[0] -= h;

	std::array<double, SIZE_NODE> xPlus2h = args;
	xPlus2h[0] += 2*h;
	std::array<double, SIZE_NODE> xMinus2h = args;
	xMinus2h[0] -= 2*h;

	std::array<double, SIZE_NODE> xPlus3h = args;
	xPlus3h[0] += 3 * h;
	std::array<double, SIZE_NODE> xMinus3h = args;
	xMinus3h[0] -= 3 * h;

	return (-function(xMinus3h) + 9 * function(xMinus2h) - 45 * function(xMinush) + 45 * function(xPlush) - 9 * function(xPlus2h) + function(xPlus3h)) / (60.0 * h);
}

inline double derivativeY(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1)
{

	std::array<double, SIZE_NODE> yPlush = args;
	yPlush[1] += h;
	std::array<double, SIZE_NODE> yMinush = args;
	yMinush[1] -= h;

	std::array<double, SIZE_NODE> yPlus2h = args;
	yPlus2h[1] += 2 * h;
	std::array<double, SIZE_NODE> yMinus2h = args;
	yMinus2h[1] -= 2 * h;

	std::array<double, SIZE_NODE> yPlus3h = args;
	yPlus3h[1] += 3 * h;
	std::array<double, SIZE_NODE> yMinus3h = args;
	yMinus3h[1] -= 3 * h;

	return (-function(yMinus3h) + 9 * function(yMinus2h) - 45 * function(yMinush) + 45 * function(yPlush) - 9 * function(yPlus2h) + function(yPlus3h)) / (60.0 * h);
}

inline double derivativeZ(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1)
{
	std::array<double, SIZE_NODE> zPlush = args;
	zPlush[2] += h;
	std::array<double, SIZE_NODE> zMinush = args;
	zMinush[2] -= h;

	std::array<double, SIZE_NODE> zPlus2h = args;
	zPlus2h[2] += 2 * h;
	std::array<double, SIZE_NODE> zMinus2h = args;
	zMinus2h[2] -= 2 * h;

	std::array<double, SIZE_NODE> zPlus3h = args;
	zPlus3h[2] += 3 * h;
	std::array<double, SIZE_NODE> zMinus3h = args;
	zMinus3h[2] -= 3 * h;

	return (-function(zMinus3h) + 9 * function(zMinus2h) - 45 * function(zMinush) + 45 * function(zPlush) - 9 * function(zPlus2h) + function(zPlus3h)) / (60.0 * h);
}

inline FiniteVector<SIZE_NODE> gradient(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args)
{
	return FiniteVector<SIZE_NODE>({ derivativeX(function, args), derivativeY(function, args), derivativeZ(function, args) });
}

inline double secondDerivativeX(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1)
{
	std::array<double, SIZE_NODE> xPlush = args;
	xPlush[0] += h;
	std::array<double, SIZE_NODE> xMinush = args;
	xMinush[0] -= h;

	std::array<double, SIZE_NODE> xPlus2h = args;
	xPlus2h[0] += 2 * h;
	std::array<double, SIZE_NODE> xMinus2h = args;
	xMinus2h[0] -= 2 * h;

	return (-2 * function(xMinus2h) + 32 * function(xMinush) - 60 * function(args) + 32 * function(xPlush) - 2 * function(xPlus2h)) / (24.0 * h * h);
}

inline double secondDerivativeY(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1)
{

	std::array<double, SIZE_NODE> yPlush = args;
	yPlush[1] += h;
	std::array<double, SIZE_NODE> yMinush = args;
	yMinush[1] -= h;

	std::array<double, SIZE_NODE> yPlus2h = args;
	yPlus2h[1] += 2 * h;
	std::array<double, SIZE_NODE> yMinus2h = args;
	yMinus2h[1] -= 2 * h;

	return (-2 * function(yMinus2h) + 32 * function(yMinush) - 60 * function(args) + 32 * function(yPlush) - 2 * function(yPlus2h)) / (24.0 * h * h);
}

inline double secondDerivativeZ(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1)
{
	std::array<double, SIZE_NODE> zPlush = args;
	zPlush[2] += h;
	std::array<double, SIZE_NODE> zMinush = args;
	zMinush[2] -= h;

	std::array<double, SIZE_NODE> zPlus2h = args;
	zPlus2h[2] += 2 * h;
	std::array<double, SIZE_NODE> zMinus2h = args;
	zMinus2h[2] -= 2 * h;

	return (-2 * function(zMinus2h) + 32 * function(zMinush) - 60 * function(args) + 32 * function(zPlush) - 2 * function(zPlus2h)) / (24.0 * h * h);
}

inline double operatorLaplace(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args)
{
	return secondDerivativeX(function, args) + secondDerivativeY(function, args) + secondDerivativeZ(function, args);
}

inline FiniteVector<SIZE_NODE> vectorProduct(const FiniteVector<SIZE_NODE>& v1, const FiniteVector<SIZE_NODE>& v2)
{
	return FiniteVector<SIZE_NODE>({ v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2] , v1[0] * v2[1] - v1[1] * v2[0] });
}

#endif



