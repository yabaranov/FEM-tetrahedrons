#pragma once

#include <vector>
#include <array>
#include <sstream>
#include <fstream>
#include <unordered_set>
#include "vector/FiniteMatrixVector.h"

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

inline std::string getFileString(const std::string& filePath)
{
	std::ifstream f;
	f.open(filePath.c_str(), std::ios::in | std::ios::binary);
	if (!f.is_open())
	{
		std::cerr << "Failed to open file: " << filePath << std::endl;
		return std::string{};
	}

	std::stringstream buffer;
	buffer << f.rdbuf();
	return buffer.str();
}

inline void removeDuplicates(std::vector<int>& v)
{
	std::unordered_set<int> s;
	auto end = std::remove_if(v.begin(), v.end(),
		[&s](int const& i) {
			return !s.insert(i).second;
		});

	v.erase(end, v.end());
}



