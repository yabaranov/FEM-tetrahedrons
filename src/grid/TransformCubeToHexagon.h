#pragma once

#include <vector>
#include <array>
#include <functional>
#include "../dimensions.h"

class TransformCubeToHexahedron
{

public:
	void Transform(std::array<double, SIZE_NODE>& node);
	TransformCubeToHexahedron(const std::array<std::pair<double, double>, SIZE_NODE>& limCube,
		const std::array<std::array<double, SIZE_NODE>, NUMBER_NODES_CUBE>& hex) :
		m_limitsCube(limCube), m_hexahedron(hex) {}

private:

	const std::array<std::array<double, SIZE_NODE>, NUMBER_NODES_CUBE>& m_hexahedron;
	const std::array<std::pair<double, double>, SIZE_NODE>& m_limitsCube;

	std::array<std::function<double(double)>, 2> m_X =
	{
		[&](double x) { return (m_limitsCube[0].second - x) / (m_limitsCube[0].second - m_limitsCube[0].first); },
		[&](double x) { return (x - m_limitsCube[0].first) / (m_limitsCube[0].second - m_limitsCube[0].first); }
	};
	std::array<std::function<double(double)>, 2> m_Y =
	{
		[&](double y) { return (m_limitsCube[1].second - y) / (m_limitsCube[1].second - m_limitsCube[1].first); },
		[&](double y) { return (y - m_limitsCube[1].first) / (m_limitsCube[1].second - m_limitsCube[1].first); }
	};
	std::array<std::function<double(double)>, 2> m_Z =
	{
		[&](double z) { return (m_limitsCube[2].second - z) / (m_limitsCube[2].second - m_limitsCube[2].first); },
		[&](double z) { return (z - m_limitsCube[2].first) / (m_limitsCube[2].second - m_limitsCube[2].first); }
	};
	double Psi(int i, const std::array<double, SIZE_NODE>& args) 
	{ return m_X[i % 2](args[0]) * m_Y[(i / 2) % 2](args[1]) * m_Z[i / 4](args[2]); }
};