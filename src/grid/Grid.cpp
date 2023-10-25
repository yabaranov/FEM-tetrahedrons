#include <fstream>
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

void Grid::TransformDomains(void)
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

void Grid::ReadGrid()
{
	std::ifstream domain("../res/input/domain.txt");
	int n = 0, m = 0;

	domain >> n >> m;
	m_supportNodes.resize(m);
	for (int i = 0; i < m; i++)
	{
		m_supportNodes[i].resize(n);
		for (int j = 0; j < n; j++)						
			domain >> m_supportNodes[i][j].x >> m_supportNodes[i][j].y;			
	}
		
	domain >> n;
	m_height.resize(n);
	for (int i = 0; i < n; i++)	
		domain >> m_height[i];
	

	domain >> n;
	m_subdomains.resize(n);

	for (int i = 0; i < n; i++)
	{
		domain >> m_subdomains[i].numberFormula;
		for (int j = 0; j < m_subdomains[i].boundaries.size(); j++)
			domain >> m_subdomains[i].boundaries[j];
	}


	std::ifstream mesh("../res/input/mesh.txt");

	m_splittingGrid.resize(SIZE_NODE);
	m_splittingGrid[0].resize(m_supportNodes[0].size() - 1);
	m_splittingGrid[1].resize(m_supportNodes.size() - 1);
	m_splittingGrid[2].resize(m_height.size() - 1);

	for (int i = 0; i < SIZE_NODE; i++)	
		for (int j = 0; j < m_splittingGrid[i].size(); j++)
			mesh >> m_splittingGrid[i][j].numIntervals >> m_splittingGrid[i][j].coefficientDischarge;

	mesh >> n;
	switch (n)
	{
		case 5:
			m_gridPattern = GRID_PATTERN::FIVE;
			break;
		case 6:
			m_gridPattern = GRID_PATTERN::SIX;
			break;
		default:
			std::cerr << "Error grid pattern!";
			break;
	}

	mesh >> n;
	switch (n)
	{
	case 0:
		m_gridNesting = NESTING::NON_NESTED;
		break;
	case 1:
		m_gridNesting = NESTING::NESTED;
		break;
	case 2:
		m_gridNesting = NESTING::DOUBLE_NESTED;
		break;
	default:
		std::cerr << "Error grid nesting!";
		break;
	}

	std::ifstream boundaryConditions("../res/input/boundary сonditions.txt");
	boundaryConditions >> n;
	m_BC.resize(n);
	for (int i = 0; i < n; i++)
	{
		boundaryConditions >> m;
		switch (m)
		{
		case 1:
			m_BC[i].typeBC = TYPE_BOUNDARY_CONDITION::FIRST;
			break;
		case 2:
			m_BC[i].typeBC = TYPE_BOUNDARY_CONDITION::SECOND;
			break;

		case 3:
			m_BC[i].typeBC = TYPE_BOUNDARY_CONDITION::THIRD;
			break;
		default:
			std::cerr << "Error type boundary condition!";
			break;
		}
		boundaryConditions >> m_BC[i].numberFormula;
		for (int j = 0; j < m_BC[i].boundaries.size(); j++)
			boundaryConditions >> m_BC[i].boundaries[j];
	}

}

void Grid::WriteGrid() const
{
	std::ofstream gridNodes("../res/output/nodes.txt");

	for (const auto& node : m_nodes)
		gridNodes << node[0] << " " << node[1] << " " << node[2] << "\n";

	std::ofstream gridElements("../res/output/elements.txt");

	for (const auto& elem : m_elements)
		gridElements << elem.vertexes[0] << " " << elem.vertexes[1] << " " <<
		elem.vertexes[2] << " " << elem.vertexes[3] << "\n";
}

int Grid::GetNumberNeighboringNodes() const
{
	return (m_gridPattern == GRID_PATTERN::FIVE) ? 
		NUMBER_NEIGHBORING_NODES::NEIGHBORS_FIVE_GRID_PATTERN : NUMBER_NEIGHBORING_NODES::NEIGHBORS_SIX_GRID_PATTERN;
}

