#include <fstream>
#include <rapidjson/document.h>
#include <rapidjson/error/en.h>
#include "Grid.h"
#include "TransformCubeToHexagon.h"

void Grid::CreateGrid()
{
	TransformDomains();
	std::array<std::vector<double>, SIZE_NODE> xyz = FormCubicGrid();
	std::vector<int> missingNodes = FormNodes(xyz);
	FormElements(xyz, missingNodes);
	//FormBC(xyz);	
}

void Grid::TransformDomains()
{
	for (int i = 0; i < m_subdomains.size(); i++)
	{
		const auto& bound = m_subdomains[i].boundaries;

		if (!QuadrilateralIsRectangle({ m_supportNodes[bound[2]][bound[0]], m_supportNodes[bound[2]][bound[1]],
			m_supportNodes[bound[3]][bound[0]], m_supportNodes[bound[3]][bound[1]] }))
		{
			std::array<std::array<double, SIZE_NODE>, NUMBER_NODES_CUBE> hexahedron;

			std::array<int, 2> ind_1 = { 2, 3 };
			std::array<int, 2> ind_2 = { 0, 1 };
			std::array<int, 2> ind_3 = { 4, 5 };

			for (int j = 0; j < hexahedron.size(); j++)
			{
				int n_1 = (j / 2) % 2;
				int n_2 = j % 2;
				int n_3 = j / 4;
				hexahedron[j][0] = m_supportNodes[bound[ind_1[n_1]]][bound[ind_2[n_2]]].x;
				hexahedron[j][1] = m_supportNodes[bound[ind_1[n_1]]][bound[ind_2[n_2]]].y;
				hexahedron[j][2] = m_height[bound[ind_3[n_3]]];
			}
			m_hexahedronSubdomains.push_back(hexahedron);
			m_numbersCurvedSubdomains.push_back(i);
		}

		m_cubicSubdomains.push_back
		({
			std::pair<double, double>{ std::min(m_supportNodes[bound[2]][bound[0]].x, m_supportNodes[bound[3]][bound[0]].x),
			std::max(m_supportNodes[bound[2]][bound[1]].x, m_supportNodes[bound[3]][bound[1]].x) },
			std::pair<double, double>{ std::min(m_supportNodes[bound[2]][bound[0]].y, m_supportNodes[bound[2]][bound[1]].y),
			std::max(m_supportNodes[bound[3]][bound[0]].y, m_supportNodes[bound[3]][bound[1]].y) },
			std::pair<double, double>{ m_height[bound[4]],
			m_height[bound[5]] } 
		});
	}
}

bool Grid::QuadrilateralIsRectangle(std::array<Point2D, NUMBER_NODES_CUBE / 2> quadrilateral) const
{
	if (quadrilateral[0].x != quadrilateral[2].x || quadrilateral[1].x != quadrilateral[3].x
		|| quadrilateral[0].y != quadrilateral[1].y || quadrilateral[2].y != quadrilateral[3].y)
		return false;

	return true;
}

std::array<std::vector<double>, SIZE_NODE> Grid::FormCubicGrid() const
{	
	//вычисление лимитов сетки
	std::array<std::vector<double>, SIZE_NODE> m_limitsСubicGrid = CalculationLimitsСubicGrid();

	// построение кубической сетки
	std::array<std::vector<double>, SIZE_NODE> xyz;

	for (int k = 0; k < xyz.size(); k++)
	{
		for (int i = 0; i < m_splittingGrid[k].size(); i++)
		{
			double b1 = 1.0;
			if (m_splittingGrid[k][i].coefficientDischarge != 1)
				b1 = (m_limitsСubicGrid[k][i + 1] - m_limitsСubicGrid[k][i]) *
				(1 - m_splittingGrid[k][i].coefficientDischarge) / (1 - pow(m_splittingGrid[k][i].coefficientDischarge, m_gridNesting * m_splittingGrid[k][i].numIntervals));
			else
				b1 = (m_limitsСubicGrid[k][i + 1] - m_limitsСubicGrid[k][i]) / (m_gridNesting * m_splittingGrid[k][i].numIntervals);

			xyz[k].push_back(m_limitsСubicGrid[k][i]);
			double degree = 1;
			for (int j = 0; j < m_gridNesting * m_splittingGrid[k][i].numIntervals - 1; j++)
			{
				xyz[k].push_back(xyz[k].back() + b1 * degree);
				degree *= m_splittingGrid[k][i].coefficientDischarge;
			}			
		}

		xyz[k].push_back(m_limitsСubicGrid[k][m_splittingGrid[k].size()]);
	}
		
	return xyz;
}

std::array<std::vector<double>, SIZE_NODE> Grid::CalculationLimitsСubicGrid(void) const
{
	std::array<std::vector<double>, SIZE_NODE> m_limitsСubicGrid;

	for (int i = 0; i < m_supportNodes[0].size(); i++)
	{
		double max = std::abs(m_supportNodes[0][i].x); int index = 0;
		for (int j = 0; j < m_supportNodes.size(); j++)
			if (std::abs(m_supportNodes[j][i].x) > max)
			{
				max = std::abs(m_supportNodes[j][i].x);
				index = j;
			}
		m_limitsСubicGrid[0].push_back(m_supportNodes[index][i].x);
	}

	for (int i = 0; i < m_supportNodes.size(); i++)
	{
		double max = std::abs(m_supportNodes[i][0].y); int index = 0;
		for (int j = 0; j < m_supportNodes[0].size(); j++)
			if (std::abs(m_supportNodes[i][j].y) > max)
			{
				max = std::abs(m_supportNodes[i][j].y);
				index = j;
			}

		m_limitsСubicGrid[1].push_back(m_supportNodes[i][index].y);
	}

	m_limitsСubicGrid[2] = m_height;

	return m_limitsСubicGrid;
}

std::pair<bool, int> Grid::InDomain(const std::array<double, SIZE_NODE>& node) const
{
	for (int i = 0; i < m_cubicSubdomains.size(); i++)
	{
		bool in = true;
		for (int j = 0; j < m_cubicSubdomains[i].size(); j++)
			if (!(m_cubicSubdomains[i][j].first <= node[j] && node[j] <= m_cubicSubdomains[i][j].second))
				in = false;
		if(in) return { true, i };
	}
	return { false, -1 };
}

std::vector<int> Grid::FormNodes(const std::array<std::vector<double>, SIZE_NODE>& xyz)
{
	std::vector<int> missingNodes;
	int numberMissingNodes = 0;

	auto&& [x, y, z] = xyz;

	for (int k = 0; k < z.size(); k++)
		for (int j = 0; j < y.size(); j++)
			for (int i = 0; i < x.size(); i++)
			{				
				std::array<double, SIZE_NODE> node = { x[i], y[j], z[k] };
				std::pair<bool, int> insideDomain = InDomain(node);
				//если узел принадлежит какой-то кубической области
				if (insideDomain.first)
				{
					//если область приведённая к кубической из шестигранной
					int index = binarySearch(m_numbersCurvedSubdomains, insideDomain.second, 0, m_numbersCurvedSubdomains.size() - 1);
					if (index >= 0)
					{
						//отображаем узел из нее в изначальную область
						TransformCubeToHexahedron isoparametricMap(m_cubicSubdomains[insideDomain.second],
							m_hexahedronSubdomains[index]);
						isoparametricMap.Transform(node);
					}
					m_nodes.push_back(node);
				}
				else
					numberMissingNodes++;
				missingNodes.push_back(numberMissingNodes);
			}
	return missingNodes;
}
				
void Grid::FormElements(const std::array<std::vector<double>, SIZE_NODE>& xyz, std::vector<int> missingNodes)
{
	auto&& [x, y, z] = xyz;

	for (int k = 0; k < z.size() - 1; k++)
	{
		int kxy_0 = k * x.size() * y.size();
		int kxy_1 = (k + 1) * x.size() * y.size();

		for (int j = 0; j < y.size() - 1; j++)
		{
			int jx_0 = j * x.size();
			int jx_1 = (j + 1) * x.size();

			for (int i = 0; i < x.size() - 1; i++)
			{
				//проверка попадания центра куба в кубические области
				std::pair<bool, int> insideDomain = InDomain(
				{ (x[i] + x[i + 1]) / 2.0, (y[j] + y[j + 1]) / 2.0, (z[k] + z[k + 1]) / 2.0 });

				//разбиение попавшего в расчётную область куба на тетраэдры
				std::array<int, NUMBER_NODES_CUBE> v =
				{
					kxy_0 + jx_0 + i,
					kxy_0 + jx_0 + i + 1,
					kxy_0 + jx_1 + i,
					kxy_0 + jx_1 + i + 1 ,
					kxy_1 + jx_0 + i,
					kxy_1 + jx_0 + i + 1,
					kxy_1 + jx_1 + i,
					kxy_1 + jx_1 + i + 1
				};

				//коррекция номеров вершин с учётом пропущенных вершин
				for (int h = 0; h < v.size(); h++)				
					v[h] -= missingNodes[v[h]];
				

				if (insideDomain.first)
				{
					auto numFormula = m_subdomains[insideDomain.second].numberFormula;
					if (m_gridPattern == GRID_PATTERN::FIVE)
					{
						if ((i + j + k) % 2 == 0)
						{
							m_elements.push_back({ { v[0], v[1], v[3], v[5] }, numFormula });
							m_elements.push_back({ { v[0], v[2], v[3], v[6] }, numFormula });
							m_elements.push_back({ { v[3], v[5], v[6], v[7] }, numFormula });
							m_elements.push_back({ { v[0], v[4], v[5], v[6] }, numFormula });
							m_elements.push_back({ { v[0], v[3], v[5], v[6] }, numFormula });
						}
						else
						{
							m_elements.push_back({ { v[0], v[1], v[2], v[4] }, numFormula });
							m_elements.push_back({ { v[1], v[2], v[3], v[7] }, numFormula });
							m_elements.push_back({ { v[2], v[4], v[6], v[7] }, numFormula });
							m_elements.push_back({ { v[1], v[4], v[5], v[7] }, numFormula });
							m_elements.push_back({ { v[1], v[2], v[4], v[7] }, numFormula });
						}
					}
					else
					{
						m_elements.push_back({{ v[0], v[1], v[2], v[6] }, numFormula });
						m_elements.push_back({{ v[0], v[1], v[4], v[6] }, numFormula });
						m_elements.push_back({{ v[1], v[4], v[5], v[6] }, numFormula });
						m_elements.push_back({{ v[1], v[2], v[3], v[6] }, numFormula });
						m_elements.push_back({{ v[1], v[3], v[6], v[7] }, numFormula });
						m_elements.push_back({{ v[1], v[5], v[6], v[7] }, numFormula });
					}
				}											
			}
		}
	}

}

//void Grid::FormBC(const std::array<std::vector<double>, SIZE_NODE>& xyz)
//{
//	auto&& [x, y, z] = xyz;
//
//	switch (m_boundaryConditions[0])
//	{
//	case 1:
//		for (int k = 0; k < z.size(); k++)
//			for (int i = 0; i < x.size(); i++)
//				m_BC_1.push_back(k * x.size() * y.size() + i);
//		break;
//	case 2:
//	case 3:
//		std::vector<Grid::Edge>&BC = (m_boundaryConditions[0] == 2) ? m_BC_2 : m_BC_3;
//		for (int k = 0; k < z.size() - 1; k++)
//		{
//			int kxy_0 = k * x.size() * y.size();
//			int kxy_1 = (k + 1) * x.size() * y.size();
//			for (int i = 0; i < x.size() - 1; i++)
//			{
//				std::array<int, 8> v =
//				{
//					kxy_0 + i,
//					kxy_0 + i + 1,
//					kxy_0 + x.size() + i,
//					kxy_0 + x.size() + i + 1 ,
//					kxy_1 + i,
//					kxy_1 + i + 1,
//					kxy_1 + x.size() + i,
//					kxy_1 + x.size() + i + 1
//				};
//				
//				if ((i + k) % 2 == 0)
//				{
//					BC.push_back({ {v[0], v[4], v[5]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + i) + 3 });
//					BC.push_back({ {v[0], v[1], v[5]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + i) + 0 });
//				}
//				else
//				{
//					BC.push_back({ {v[0], v[1], v[4]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + i) + 0 });
//					BC.push_back({ {v[1], v[4], v[5]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + i) + 3 });
//				}				
//			}
//		}
//		break;
//	}
//
//	switch (m_boundaryConditions[1])
//	{
//	case 1:
//		for (int k = 0; k < z.size(); k++)
//			for (int i = 0; i < x.size(); i++)
//				m_BC_1.push_back(k * x.size() * y.size() + (y.size() - 1) * x.size() + i);
//		break;
//	case 2:
//	case 3:
//		std::vector<Grid::Edge>&BC = (m_boundaryConditions[1] == 2) ? m_BC_2 : m_BC_3;
//		for (int k = 0; k < z.size() - 1; k++)
//		{
//			int kxy_0 = k * x.size() * y.size();
//			int kxy_1 = (k + 1) * x.size() * y.size();
//
//			for (int i = 0; i < x.size() - 1; i++)
//			{
//				std::array<int, 8> v =
//				{
//					kxy_0 + (y.size() - 2) * x.size() + i,
//					kxy_0 + (y.size() - 2) * x.size() + i + 1,
//					kxy_0 + (y.size() - 1) * x.size() + i,
//					kxy_0 + (y.size() - 1) * x.size() + i + 1 ,
//					kxy_1 + (y.size() - 2) * x.size() + i,
//					kxy_1 + (y.size() - 2) * x.size() + i + 1,
//					kxy_1 + (y.size() - 1) * x.size() + i,
//					kxy_1 + (y.size() - 1) * x.size() + i + 1
//				};
//				
//				if ((i + y.size() - 2 + k) % 2 == 0)
//				{
//					BC.push_back({ {v[2], v[3], v[6]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + (y.size() - 2) * (x.size() - 1) + i) + 1 });
//					BC.push_back({ {v[3], v[6], v[7]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + (y.size() - 2) * (x.size() - 1) + i) + 2 });
//				}
//				else
//				{
//					BC.push_back({ {v[2], v[3], v[7]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + (y.size() - 2) * (x.size() - 1) + i) + 1 });
//					BC.push_back({ {v[2], v[6], v[7]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + (y.size() - 2) * (x.size() - 1) + i) + 2 });
//				}
//				
//			}
//		}
//
//		break;
//	}
//
//	switch (m_boundaryConditions[2])
//	{
//	case 1:
//		for (int k = 0; k < z.size(); k++)
//			for (int j = 0; j < y.size(); j++)
//				m_BC_1.push_back(k * x.size() * y.size() + j * x.size());
//		break;
//	case 2:
//	case 3:
//		std::vector<Grid::Edge>&BC = (m_boundaryConditions[2] == 2) ? m_BC_2 : m_BC_3;
//		for (int k = 0; k < z.size() - 1; k++)
//		{
//			int kxy_0 = k * x.size() * y.size();
//			int kxy_1 = (k + 1) * x.size() * y.size();
//
//			for (int j = 0; j < y.size() - 1; j++)
//			{
//				int jx_0 = j * x.size();
//				int jx_1 = (j + 1) * x.size();
//				std::array<int, 8> v =
//				{
//					kxy_0 + jx_0,
//					kxy_0 + jx_0 + 1,
//					kxy_0 + jx_1,
//					kxy_0 + jx_1 + 1 ,
//					kxy_1 + jx_0,
//					kxy_1 + jx_0 + 1,
//					kxy_1 + jx_1,
//					kxy_1 + jx_1 + 1
//				};
//				
//				if ((j + k) % 2 == 0)
//				{
//					BC.push_back({ {v[0], v[4], v[6]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1)) + 3 });
//					BC.push_back({ {v[0], v[2], v[6]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1)) + 1 });
//				}
//				else
//				{
//					BC.push_back({ {v[0], v[2], v[4]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1)) + 0 });
//					BC.push_back({ {v[2], v[4], v[6]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1)) + 2 });
//				}
//				
//			}
//		}
//		break;
//	}
//
//	switch (m_boundaryConditions[3])
//	{
//	case 1:
//		for (int k = 0; k < z.size(); k++)
//			for (int j = 0; j < y.size(); j++)
//				m_BC_1.push_back(k * x.size() * y.size() + j * x.size() + x.size() - 1);
//		break;
//	case 2:
//	case 3:
//		std::vector<Grid::Edge>&BC = (m_boundaryConditions[3] == 2) ? m_BC_2 : m_BC_3;
//		for (int k = 0; k < z.size() - 1; k++)
//		{
//			int kxy_0 = k * x.size() * y.size();
//			int kxy_1 = (k + 1) * x.size() * y.size();
//
//			for (int j = 0; j < y.size() - 1; j++)
//			{
//				int jx_0 = j * x.size();
//				int jx_1 = (j + 1) * x.size();
//
//				std::array<int, 8> v =
//				{
//					kxy_0 + jx_0 + x.size() - 2,
//					kxy_0 + jx_0 + x.size() - 1,
//					kxy_0 + jx_1 + x.size() - 2,
//					kxy_0 + jx_1 + x.size() - 1,
//					kxy_1 + jx_0 + x.size() - 2,
//					kxy_1 + jx_0 + x.size() - 1,
//					kxy_1 + jx_1 + x.size() - 2,
//					kxy_1 + jx_1 + x.size() - 1
//				};
//				
//				if ((x.size() - 2 + j + k) % 2 == 0)
//				{
//					BC.push_back({ {v[1], v[3], v[5]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + x.size() - 2) + 0 });
//					BC.push_back({ {v[3], v[5], v[7]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + x.size() - 2) + 2 });
//				}
//				else
//				{
//					BC.push_back({ {v[1], v[5], v[7]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + x.size() - 2) + 3 });
//					BC.push_back({ {v[1], v[3], v[7]}, 5 * (k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + x.size() - 2) + 1 });
//				}
//				
//			}
//		}
//		break;
//	}
//
//	switch (m_boundaryConditions[4])
//	{
//	case 1:
//		for (int j = 0; j < y.size(); j++)
//			for (int i = 0; i < x.size(); i++)
//				m_BC_1.push_back(j * x.size() + i);
//		break;
//	case 2:
//	case 3:
//		std::vector<Grid::Edge>&BC = (m_boundaryConditions[4] == 2) ? m_BC_2 : m_BC_3;
//		for (int j = 0; j < y.size() - 1; j++)
//		{
//			int jx_0 = j * x.size();
//			int jx_1 = (j + 1) * x.size();
//
//			for (int i = 0; i < x.size() - 1; i++)
//			{
//				std::array<int, 8> v =
//				{
//					jx_0 + i,
//					jx_0 + i + 1,
//					jx_1 + i,
//					jx_1 + i + 1 ,
//					x.size() * y.size() + jx_0 + i,
//					x.size() * y.size() + jx_0 + i + 1,
//					x.size() * y.size() + jx_1 + i,
//					x.size() * y.size() + jx_1 + i + 1
//				};
//				
//				if ((i + j) % 2 == 0)
//				{
//					BC.push_back({ {v[0], v[1], v[3]}, 5 * (j * (x.size() - 1) + i) + 0 });
//					BC.push_back({ {v[0], v[2], v[3]}, 5 * (j * (x.size() - 1) + i) + 1 });
//				}
//				else
//				{
//					BC.push_back({ {v[0], v[1], v[2]}, 5 * (j * (x.size() - 1) + i) + 0 });
//					BC.push_back({ {v[1], v[2], v[3]}, 5 * (j * (x.size() - 1) + i) + 1 });
//				}
//				
//			}
//		}
//		break;
//	}
//
//	switch (m_boundaryConditions[5])
//	{
//	case 1:
//		for (int j = 0; j < y.size(); j++)
//			for (int i = 0; i < x.size(); i++)
//				m_BC_1.push_back((z.size() - 1) * x.size() * y.size() + j * x.size() + i);
//		break;
//	case 2:
//	case 3:
//		std::vector<Grid::Edge>&BC = (m_boundaryConditions[5] == 2) ? m_BC_2 : m_BC_3;
//		for (int j = 0; j < y.size() - 1; j++)
//		{
//			int jx_0 = j * x.size();
//			int jx_1 = (j + 1) * x.size();
//
//			for (int i = 0; i < x.size() - 1; i++)
//			{
//				std::array<int, 8> v =
//				{
//					(z.size() - 2) * x.size() * y.size() + jx_0 + i,
//					(z.size() - 2) * x.size() * y.size() + jx_0 + i + 1,
//					(z.size() - 2) * x.size() * y.size() + jx_1 + i,
//					(z.size() - 2) * x.size() * y.size() + jx_1 + i + 1 ,
//					(z.size() - 1) * x.size() * y.size() + jx_0 + i,
//					(z.size() - 1) * x.size() * y.size() + jx_0 + i + 1,
//					(z.size() - 1) * x.size() * y.size() + jx_1 + i,
//					(z.size() - 1) * x.size() * y.size() + jx_1 + i + 1
//				};
//				
//				if ((i + j + z.size() - 2) % 2 == 0)
//				{
//					BC.push_back({ {v[4], v[5], v[6]}, 5 * ((z.size() - 2) * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + i) + 3 });
//					BC.push_back({ {v[5], v[6], v[7]}, 5 * ((z.size() - 2) * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + i) + 2 });
//				}
//				else
//				{
//					BC.push_back({ {v[4], v[5], v[7]}, 5 * ((z.size() - 2) * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + i) + 3 });
//					BC.push_back({ {v[4], v[6], v[7]}, 5 * ((z.size() - 2) * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + i) + 2 });
//				}
//				
//			}
//		}
//		break;
//	}
//}

bool Grid::ReadGridJSON()
{
	std::string JSONPath = "../res/input/grid.json";
	const std::string JSONString = getFileString(JSONPath);

	if (JSONString.empty())
	{
		std::cerr << "No JSON resources file!" << std::endl;
		return false;
	}

	rapidjson::Document document;
	rapidjson::ParseResult parseResult = document.Parse(JSONString.c_str());
	if (!parseResult)
	{
		std::cerr << "JSON parse error: " << rapidjson::GetParseError_En(parseResult.Code()) << "(" << parseResult.Offset() << ")" << std::endl;
		std::cerr << "In JSON file: " << JSONPath << std::endl;
		return false;
	}

	auto domainIt = document.FindMember("domain");
	if (domainIt != document.MemberEnd())
	{
		const auto supportNodesArray = domainIt->value["supportNodes"].GetArray();
		m_supportNodes.resize(supportNodesArray.Size());
		for (size_t i = 0; i < supportNodesArray.Size(); i++)
		{
			m_supportNodes[i].resize(supportNodesArray[i].Size());
			for (size_t j = 0; j < supportNodesArray[i].Size(); j++)
			{
				const auto x = supportNodesArray[i][j]["x"].GetDouble();
				const auto y = supportNodesArray[i][j]["y"].GetDouble();
				m_supportNodes[i][j] = Point2D{ x, y };
			}
		}

		const auto heightArray = domainIt->value["height"].GetArray();
		m_height.resize(heightArray.Size());
		for (size_t i = 0; i < heightArray.Size(); i++)
			m_height[i] = heightArray[i].GetDouble();

		const auto subdomainsArray = domainIt->value["subdomains"].GetArray();
		m_subdomains.resize(subdomainsArray.Size());
		for (size_t i = 0; i < subdomainsArray.Size(); i++)
		{
			const auto numberFormula = subdomainsArray[i]["numberFormula"].GetInt();

			const auto subdomainsArrayBoundaries = subdomainsArray[i]["boundaries"].GetArray();
			std::array<int, SIZE_SUBDOMAIN> boundaries;
			for (size_t j = 0; j < subdomainsArrayBoundaries.Size(); j++)
				boundaries[j] = subdomainsArrayBoundaries[j].GetUint();

			std::array<int, SIZE_SUBDOMAIN / 2> rlimits{ m_supportNodes[0].size(), m_supportNodes.size(), m_height.size() };
			for (size_t j = 0; j < boundaries.size(); j++)
				if (boundaries[j] < 0 || boundaries[j] >= rlimits[j / 2])
				{
					std::cerr << "Number boundaries in subdomains[" + std::to_string(i) + "] is wrong!" << std::endl;
					return false;
				}

			m_subdomains[i] = Subdomain{ numberFormula,  boundaries };
		}
	}

	auto meshIt = document.FindMember("mesh");
	if (meshIt != document.MemberEnd())
	{
		const auto splittingGridArray = meshIt->value["splittingGrid"].GetArray();
		for (size_t i = 0; i < splittingGridArray.Size(); i++)
		{
			m_splittingGrid[i].resize(splittingGridArray[i].Size());
			for (size_t j = 0; j < splittingGridArray[i].Size(); j++)
			{
				auto numIntervals = splittingGridArray[i][j]["numIntervals"].GetInt();

				auto coefficientDischarge = splittingGridArray[i][j]["coefficientDischarge"].GetDouble();
				if (coefficientDischarge < 0)
				{
					std::cerr << "Coefficient discharge in splittingGrid[" + std::to_string(i) + "] < 0!" << std::endl;
					return false;
				}
				m_splittingGrid[i][j] = Splitting{ numIntervals, coefficientDischarge };
			}
		}

		const auto gridPattern = meshIt->value["gridPattern"].GetInt();
		switch (gridPattern)
		{
			case 5:
				m_gridPattern = GRID_PATTERN::FIVE; break;
			case 6:
				m_gridPattern = GRID_PATTERN::SIX; break;
			default:
				std::cerr << "Error grid pattern!"; return false;
		}

		const auto gridNesting = meshIt->value["gridNesting"].GetInt();
		switch (gridNesting)
		{
			case 0:
				m_gridNesting = NESTING::NON_NESTED; break;
			case 1:
				m_gridNesting = NESTING::NESTED; break;
			case 2:
				m_gridNesting = NESTING::DOUBLE_NESTED; break;
			default:
				std::cerr << "Error grid nesting!"; return false;
		}
	} 

	auto boundaryConditionsIt = document.FindMember("boundaryConditions");
	if (meshIt != document.MemberEnd())
	{
		const auto boundaryConditionsArray = boundaryConditionsIt->value.GetArray();
		m_BC.resize(boundaryConditionsArray.Size());

		for (size_t i = 0; i < boundaryConditionsArray.Size(); i++)
		{
			auto typeboundaryCondition = boundaryConditionsArray[i]["typeBC"].GetInt();
			TYPE_BOUNDARY_CONDITION typeBC;

			switch (typeboundaryCondition)
			{
			case 1:
				typeBC = TYPE_BOUNDARY_CONDITION::FIRST; break;
			case 2:
				typeBC = TYPE_BOUNDARY_CONDITION::SECOND; break;
			case 3:
				typeBC = TYPE_BOUNDARY_CONDITION::THIRD; break;
			default:
				std::cerr << "Error type boundary condition in boundaryConditions[" + std::to_string(i) + "]!"; return false;
			}

			int numberFormula = boundaryConditionsArray[i]["numberFormula"].GetUint();
			const auto boundaryConditionsArrayBoundaries = boundaryConditionsArray[i]["boundaries"].GetArray();
			std::array<int, SIZE_SUBDOMAIN> boundaries;
			for (size_t j = 0; j < boundaryConditionsArrayBoundaries.Size(); j++)
				boundaries[j] = boundaryConditionsArrayBoundaries[j].GetUint();

			std::array<int, SIZE_SUBDOMAIN / 2> rlimits{ m_supportNodes[0].size(), m_supportNodes.size(), m_height.size() };
			for (size_t j = 0; j < boundaries.size(); j++)
				if (boundaries[j] < 0 || boundaries[j] >= rlimits[j / 2])
				{
					std::cerr << "Number boundaries in subdomains[" + std::to_string(i) + "] is wrong!" << std::endl;
					return false;
				}

			m_BC[i] = BoundaryCondition{ typeBC,  numberFormula, boundaries };
		}
	}
	return true;
}

bool Grid::WriteGrid() const
{
	std::string filePath = "../res/output/nodes.txt";
	std::ofstream gridNodes;
	gridNodes.open(filePath.c_str(), std::ios::out | std::ios::binary);
	if (!gridNodes.is_open())
	{
		std::cerr << "Failed to open file: " << filePath << std::endl;
		return false;
	}
	for (const auto& node : m_nodes)
		gridNodes << node[0] << " " << node[1] << " " << node[2] << "\n";


	filePath = "../res/output/elements.txt";
	std::ofstream gridElements;
	gridElements.open(filePath.c_str(), std::ios::out | std::ios::binary);
	if (!gridElements.is_open())
	{
		std::cerr << "Failed to open file: " << filePath << std::endl;
		return false;
	}
	for (const auto& elem : m_elements)
		gridElements << elem.vertexes[0] << " " << elem.vertexes[1] << " " <<
		elem.vertexes[2] << " " << elem.vertexes[3] << "\n";
}

int Grid::GetNumberNeighboringNodes() const
{
	return (m_gridPattern == GRID_PATTERN::FIVE) ? 
		NUMBER_NEIGHBORING_NODES::NEIGHBORS_FIVE_GRID_PATTERN : NUMBER_NEIGHBORING_NODES::NEIGHBORS_SIX_GRID_PATTERN;
}

