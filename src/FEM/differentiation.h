#pragma once
#include <functional>
#include "../vector/FiniteMatrixVector.h"
#include "../dimensions.h"

double derivativeX(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1);
double derivativeY(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1);
double derivativeZ(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1);
FiniteVector<SIZE_NODE> gradient(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args);

double secondDerivativeX(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1);
double secondDerivativeY(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1);
double secondDerivativeZ(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args, double h = 1);
double operatorLaplace(std::function<double(const std::array<double, SIZE_NODE>&)> function, const std::array<double, SIZE_NODE>& args);

FiniteVector<SIZE_NODE> vectorProduct(const FiniteVector<SIZE_NODE>& v1, const FiniteVector<SIZE_NODE>& v2);
