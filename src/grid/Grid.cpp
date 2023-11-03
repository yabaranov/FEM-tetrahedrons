#include <fstream>
#include <rapidjson/document.h>
#include <rapidjson/error/en.h>
#include "Grid.h"
#include "TransformCubeToHexagon.h"
#include "../functions.h"

void Grid::CreateGrid()
{
	//коррекция интервалов разбиения с учетом вложенности сетки
	for (int k = 0; k < m_splittingGrid.size(); k++)
		for (int i = 0; i < m_splittingGrid[k].size(); i++)
			m_splittingGrid[k][i].numIntervals *= m_gridNesting;

	TransformDomains();
	std::array<std::vector<double>, SIZE_NODE> xyz = FormCubicGrid();
	std::vector<int> missingNodes = FormNodes(xyz);
	std::vector<int> missingElements = FormElements(xyz, missingNodes);
	FormBC(xyz, missingNodes, missingElements);
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
	std::array<std::vector<double>, SIZE_NODE> limitsСubicGrid = CalculationLimitsСubicGrid();

	// построение кубической сетки
	std::array<std::vector<double>, SIZE_NODE> xyz;

	for (int k = 0; k < xyz.size(); k++)
	{
		for (int i = 0; i < m_splittingGrid[k].size(); i++)
		{
			double b1 = 1.0;
			if (m_splittingGrid[k][i].coefficientDischarge != 1)
				b1 = (limitsСubicGrid[k][i + 1] - limitsСubicGrid[k][i]) *
				(1 - m_splittingGrid[k][i].coefficientDischarge) / (1 - pow(m_splittingGrid[k][i].coefficientDischarge, m_splittingGrid[k][i].numIntervals));
			else
				b1 = (limitsСubicGrid[k][i + 1] - limitsСubicGrid[k][i]) / (m_splittingGrid[k][i].numIntervals);

			xyz[k].push_back(limitsСubicGrid[k][i]);
			double degree = 1;
			for (int j = 0; j < m_splittingGrid[k][i].numIntervals - 1; j++)
			{
				xyz[k].push_back(xyz[k].back() + b1 * degree);
				degree *= m_splittingGrid[k][i].coefficientDischarge;
			}			
		}

		xyz[k].push_back(limitsСubicGrid[k][m_splittingGrid[k].size()]);
	}
		
	return xyz;
}

std::array<std::vector<double>, SIZE_NODE> Grid::CalculationLimitsСubicGrid() const
{
	std::array<std::vector<double>, SIZE_NODE> limitsСubicGrid;

	for (int i = 0; i < m_supportNodes[0].size(); i++)
	{
		double max = std::abs(m_supportNodes[0][i].x); int index = 0;
		for (int j = 0; j < m_supportNodes.size(); j++)
			if (std::abs(m_supportNodes[j][i].x) > max)
			{
				max = std::abs(m_supportNodes[j][i].x);
				index = j;
			}
		limitsСubicGrid[0].push_back(m_supportNodes[index][i].x);
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

		limitsСubicGrid[1].push_back(m_supportNodes[i][index].y);
	}

	limitsСubicGrid[2] = m_height;

	return limitsСubicGrid;
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
				
std::vector<int> Grid::FormElements(const std::array<std::vector<double>, SIZE_NODE>& xyz, std::vector<int> missingNodes)
{
	std::vector<int> missingElements;
	int numberMissingElements = 0;

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
				std::array<int, NUMBER_NODES_CUBE> vCube =
				{
					kxy_0 + jx_0 + i,
					kxy_0 + jx_0 + i + 1,
					kxy_0 + jx_1 + i,
					kxy_0 + jx_1 + i + 1,
					kxy_1 + jx_0 + i,
					kxy_1 + jx_0 + i + 1,
					kxy_1 + jx_1 + i,
					kxy_1 + jx_1 + i + 1
				};

				//коррекция номеров вершин с учётом пропущенных вершин
				for (int h = 0; h < vCube.size(); h++)
					vCube[h] -= missingNodes[vCube[h]];
				
				if (insideDomain.first)
				{
					auto numFormula = m_subdomains[insideDomain.second].numberFormula;
					if (m_gridPattern == GRID_PATTERN::FIVE)
					{
						if ((i + j + k) % 2 == 0)
						{
							m_elements.push_back({ { vCube[0], vCube[1], vCube[3], vCube[5] }, numFormula });
							m_elements.push_back({ { vCube[0], vCube[2], vCube[3], vCube[6] }, numFormula });
							m_elements.push_back({ { vCube[3], vCube[5], vCube[6], vCube[7] }, numFormula });
							m_elements.push_back({ { vCube[0], vCube[4], vCube[5], vCube[6] }, numFormula });
							m_elements.push_back({ { vCube[0], vCube[3], vCube[5], vCube[6] }, numFormula });
						}
						else
						{
							m_elements.push_back({ { vCube[0], vCube[1], vCube[2], vCube[4] }, numFormula });
							m_elements.push_back({ { vCube[1], vCube[2], vCube[3], vCube[7] }, numFormula });
							m_elements.push_back({ { vCube[2], vCube[4], vCube[6], vCube[7] }, numFormula });
							m_elements.push_back({ { vCube[1], vCube[4], vCube[5], vCube[7] }, numFormula });
							m_elements.push_back({ { vCube[1], vCube[2], vCube[4], vCube[7] }, numFormula });
						}
					}
					else
					{
						m_elements.push_back({{ vCube[0], vCube[1], vCube[2], vCube[6] }, numFormula });
						m_elements.push_back({{ vCube[0], vCube[1], vCube[4], vCube[6] }, numFormula });
						m_elements.push_back({{ vCube[1], vCube[4], vCube[5], vCube[6] }, numFormula });
						m_elements.push_back({{ vCube[1], vCube[2], vCube[3], vCube[6] }, numFormula });
						m_elements.push_back({{ vCube[1], vCube[3], vCube[6], vCube[7] }, numFormula });
						m_elements.push_back({{ vCube[1], vCube[5], vCube[6], vCube[7] }, numFormula });
					}
				}	
				else
					numberMissingElements++;
				missingElements.push_back(numberMissingElements);
			}
		}
	}

	return missingElements;

}

std::array<std::pair<int, int>, SIZE_NODE> Grid::CalculationLimitsEdge(const std::array<int, SIZE_SUBDOMAIN>& boundaries)
{
	std::array<std::pair<int, int>, SIZE_NODE> limitsEdge;

	for (int i = 0; i < limitsEdge.size(); i++)
	{
		std::pair<int, int> limit = { 0, 0 };
		for (int j = 0; j < boundaries[2 * i]; j++)
			limit.first += m_splittingGrid[i][j].numIntervals;

		for (int j = 0; j < boundaries[2 * i + 1]; j++)
			limit.second += m_splittingGrid[i][j].numIntervals;

		limitsEdge[i] = limit;
	}

	return limitsEdge;
}

void Grid::FormBC(const std::array<std::vector<double>, SIZE_NODE>& xyz, const std::vector<int>& missingNodes, const std::vector<int>& missingElements)
{
	auto&& [x, y, z] = xyz;

	for (int h = 0; h < m_BC.size(); h++)
	{
		std::array<std::pair<int, int>, SIZE_NODE> limitsEdge = CalculationLimitsEdge(m_BC[h].boundaries);

		switch (m_BC[h].typeBC)
		{
			case TYPE_BOUNDARY_CONDITION::FIRST:
				
				for (int k = limitsEdge[2].first; k <= limitsEdge[2].second; k++)
					for (int j = limitsEdge[1].first; j <= limitsEdge[1].second; j++)
						for (int i = limitsEdge[0].first; i <= limitsEdge[0].second; i++)
						{
								int number = k * x.size() * y.size() + j * x.size() + i;
							//коррекция номеров вершин с учётом пропущенных вершин
							number -= missingNodes[number];
							m_BC_1.push_back(number);
						}	
				//удаление дубликатов из массива первых краевых(дубликаты расположены на рёбрах, где грани соприкасаются)
				removeDuplicates(m_BC_1);
				break;

			case TYPE_BOUNDARY_CONDITION::SECOND:
			case TYPE_BOUNDARY_CONDITION::THIRD:
				std::vector<Edge>& BC = (m_BC[h].typeBC == TYPE_BOUNDARY_CONDITION::SECOND) ? m_BC_2 : m_BC_3;

				bool reducelimit = false;
				int coordinateMatching = 0;
				//коррекция лимитов циклов
				for (int i = 0; i < limitsEdge.size(); i++)
					if (limitsEdge[i].first != limitsEdge[i].second)					
						limitsEdge[i].second--;				
					else
					{
						coordinateMatching = i;		
						if (0 < limitsEdge[i].first)
						{
							std::array<double, SIZE_NODE> leftNode{0, 0, 0};
							leftNode[i]--;
							for (int j = 0; j < leftNode.size(); j++)
								leftNode[j] = (xyz[j][limitsEdge[j].first + leftNode[j]] + xyz[j][limitsEdge[j].second + leftNode[j]]) / 2.0;
							if (InDomain(leftNode).first)
								reducelimit = true;
						}					
					}

				for (int k = limitsEdge[2].first; k <= limitsEdge[2].second; k++)
				{
					int kxy_0 = k * x.size() * y.size();
					int kxy_1 = (k + 1) * x.size() * y.size();

					for (int j = limitsEdge[1].first; j <= limitsEdge[1].second; j++)
					{
						int jx_0 = j * x.size();
						int jx_1 = (j + 1) * x.size();

						for (int i = limitsEdge[0].first; i <= limitsEdge[0].second; i++)
						{
							std::array<int, NUMBER_NODES_CUBE / 2> vRect;
							int elementOffset = k * (x.size() - 1) * (y.size() - 1) + j * (x.size() - 1) + i;
							//вычисление вершин прямоугольника, который будет разбит на грани(треугольники)
							switch (coordinateMatching)
							{
							case 0:
								vRect = { kxy_0 + jx_0 + i, kxy_0 + jx_1 + i, kxy_1 + jx_0 + i, kxy_1 + jx_1 + i };
								if(reducelimit) elementOffset--;								
								break;
							case 1:
								vRect = { kxy_0 + jx_0 + i, kxy_0 + jx_0 + i + 1, kxy_1 + jx_0 + i, kxy_1 + jx_0 + i + 1 };
								if (reducelimit) elementOffset -= x.size() - 1;								
								break;
							case 2:
								vRect = { kxy_0 + jx_0 + i, kxy_0 + jx_0 + i + 1, kxy_0 + jx_1 + i, kxy_0 + jx_1 + i + 1 };
								if (reducelimit) elementOffset -= (x.size() - 1) * (y.size() - 1);
								break;
							}

							//коррекция номеров вершин с учётом пропущенных вершин
							for (int l = 0; l < vRect.size(); l++)
								vRect[l] -= missingNodes[vRect[l]];
							
							//коррекция номера элемента-куба с учётом пропущенных элементов-кубов
							elementOffset -= missingElements[elementOffset];
							//каждый куб был разбит на тетраэдры
							elementOffset *= m_gridPattern;
						
							std::array<int, SIZE_EDGE> vertexes;
							if (m_gridPattern == GRID_PATTERN::FIVE)
							{
								if ((i + j + k) % 2 == 0)
								{
									vertexes = { vRect[0], vRect[1], vRect[3] };
									BC.push_back({ vertexes, SearchElement(vertexes, elementOffset, elementOffset + m_gridPattern - 1) });
									vertexes = { vRect[0], vRect[2], vRect[3] };
									BC.push_back({ vertexes, SearchElement(vertexes, elementOffset, elementOffset + m_gridPattern - 1) });
								}							
								else
								{
									vertexes = { vRect[0], vRect[1], vRect[2] };
									BC.push_back({ vertexes, SearchElement(vertexes, elementOffset, elementOffset + m_gridPattern - 1) });
									vertexes = { vRect[1], vRect[2], vRect[3] };
									BC.push_back({ vertexes, SearchElement(vertexes, elementOffset, elementOffset + m_gridPattern - 1) });
								}
							}
							else
							{
								vertexes = { vRect[0], vRect[1], vRect[2] };
								BC.push_back({ vertexes, SearchElement(vertexes, elementOffset, elementOffset + m_gridPattern - 1) });
								vertexes = { vRect[1], vRect[2], vRect[3] };
								BC.push_back({ vertexes, SearchElement(vertexes, elementOffset, elementOffset + m_gridPattern - 1) });
							}

						}
					}
				}
												
				break;
		}
	}
	
}

int Grid::SearchElement(const std::array<int, SIZE_EDGE>& vertexes, int l, int r)
{
	for (int i = l; i <= r; i++)
	{
		bool edgeInElement = true;
		for (int j = 0; j < vertexes.size(); j++)		
			if (binarySearch(m_elements[i].vertexes, vertexes[j], 0, m_elements[i].vertexes.size() - 1) == -1)
				edgeInElement = false;		
		if (edgeInElement) return i;
	}
	return -1;
}

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
				boundaries[j] = subdomainsArrayBoundaries[j].GetInt();

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

			const auto boundaryConditionsArrayBoundaries = boundaryConditionsArray[i]["boundaries"].GetArray();
			std::array<int, SIZE_SUBDOMAIN> boundaries;
			for (size_t j = 0; j < boundaryConditionsArrayBoundaries.Size(); j++)
				boundaries[j] = boundaryConditionsArrayBoundaries[j].GetInt();

			std::array<int, SIZE_SUBDOMAIN / 2> rlimits{ m_supportNodes[0].size(), m_supportNodes.size(), m_height.size() };
			for (size_t j = 0; j < boundaries.size(); j++)
				if (boundaries[j] < 0 || boundaries[j] >= rlimits[j / 2])
				{
					std::cerr << "Number boundaries in subdomains[" + std::to_string(i) + "] is wrong!" << std::endl;
					return false;
				}

			m_BC[i] = BoundaryCondition{ typeBC, boundaries };
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

