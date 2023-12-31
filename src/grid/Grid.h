#pragma once

#include <vector>
#include <array>
#include <functional>
#include "../dimensions.h"

class Grid
{
public:
	Grid() {}
	~Grid() {}
	bool ReadGridJSON();
	void CreateGrid();
	bool WriteGrid() const;
	
	enum TYPE_BOUNDARY_CONDITION { FIRST = 1, SECOND = 2, THIRD = 3 };

	const std::array<double, SIZE_NODE>& GetNode(int i) const { return m_nodes[i]; }

	struct Element
	{
		std::array<int, SIZE_ELEMENT> vertexes;
		int numberFormula;
	};
	const std::vector<Element>& Elements() const { return m_elements; }
	const Element& GetElement(int i) const { return m_elements[i]; }

	struct BoundaryEdge
	{		
		std::array<int, SIZE_EDGE> vertexes;
		int elemNum;
	};
	const BoundaryEdge& GetBoundaryEdgeBC_2(int i) const { return m_BC_2[i]; }
	const BoundaryEdge& GetBoundaryEdgeBC_3(int i) const { return m_BC_3[i]; }
	int GetNumberNodeBC_1(int i) const { return m_BC_1[i]; }

	struct Edge
	{
		std::array<int, SIZE_EDGE> vertexes;
		std::array<int, 2> elemNums;
		bool operator!=(const Edge& edge) const { return edge.vertexes != vertexes; }
		bool operator==(const Edge& edge) const { return !(edge != *this); }
	};
	const Edge& GetEdge(int i) const { return m_edges[i]; }
	int GetNumberEdge(const Edge& edge) const;
	const std::array<int, NUMBER_EDGES_ELEMENT>& GetNumbersEdges(int numElem) const;

	struct Subdomain
	{
		int numberFormula;
		std::array<int, SIZE_SUBDOMAIN> boundaries;
	};
	Subdomain GetSubdomain(int i) { return m_subdomains[i]; }


	int SizeElements() const { return m_elements.size(); }
	int SizeNodes() const { return m_nodes.size(); }
	int SizeBoundaryEdgesBC_2() const { return m_BC_2.size(); }
	int SizeBoundaryEdgesBC_3() const { return m_BC_3.size(); }
	int SizeNodesBC_1() const { return m_BC_1.size(); }

	int GetNumberNeighboringNodes() const;

	enum NESTING { NON_NESTED = 1, NESTED = 2, DOUBLE_NESTED = 4 };
	NESTING GetGridNesting() const { return m_gridNesting; }

private:
	enum NUMBER_NEIGHBORING_NODES { NEIGHBORS_FIVE_GRID_PATTERN = 18, NEIGHBORS_SIX_GRID_PATTERN = 14 };

	enum GRID_PATTERN { FIVE = 5, SIX = 6 };
	GRID_PATTERN m_gridPattern = GRID_PATTERN::FIVE;

	
	NESTING m_gridNesting = NESTING::NON_NESTED;

	std::vector<std::array<double, SIZE_NODE>> m_nodes;
	std::vector<Element> m_elements;
	std::vector<BoundaryEdge> m_BC_2;
	std::vector<BoundaryEdge> m_BC_3;
	std::vector<int> m_BC_1;	
	std::vector<Edge> m_edges;
	std::vector<Subdomain> m_subdomains;

	struct Splitting
	{
		int numIntervals;
		double coefficientDischarge;		
	};

	std::array<std::vector<Splitting>, SIZE_NODE> m_splittingGrid;

	struct Point2D
	{
		double x;
		double y;
	};
	std::vector<std::vector<Point2D>> m_supportNodes;
	std::vector<double> m_height;

	struct BoundaryCondition
	{
		TYPE_BOUNDARY_CONDITION typeBC;
		std::array<int, SIZE_SUBDOMAIN> boundaries;		
	};

	std::vector<BoundaryCondition> m_BC;

	std::vector<std::array<std::array<double, SIZE_NODE>, NUMBER_NODES_CUBE>> m_trulyHexahedronSubdomains;
	std::vector<std::array<std::pair<double, double>, SIZE_NODE>> m_cubicSubdomains;
	std::vector<int> m_numbersCurvedSubdomains;

	std::vector<std::array<std::array<double, SIZE_NODE>, NUMBER_NODES_CUBE>> m_hexahedronSubdomains;

	void TransformDomains();
	bool QuadrilateralIsRectangle(std::array<Point2D, NUMBER_NODES_CUBE / 2> quadrilateral) const;
	std::array<std::vector<double>, SIZE_NODE> FormCubicGrid() const;
	std::array<std::vector<double>, SIZE_NODE> CalculationLimits�ubicGrid() const;
	std::pair<bool, int> InDomain(const std::array<double, SIZE_NODE>& node) const;

	std::vector<int> FormNodes(const std::array<std::vector<double>, SIZE_NODE>& xyz);
	std::vector<int> FormElements(const std::array<std::vector<double>, SIZE_NODE>& xyz, std::vector<int> missingNodes);
	void FormBC(const std::array<std::vector<double>, SIZE_NODE>& xyz, const std::vector<int>& missingNodes, const std::vector<int>& missingElements);
	std::array<std::pair<int, int>, SIZE_NODE> CalculationLimitsBoundaryEdge(const std::array<int, SIZE_SUBDOMAIN>& boundaries);
	int SearchElement(const std::array<int, SIZE_EDGE>& vertexes, int l, int r);
	void FormEdges();

};