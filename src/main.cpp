#include "grid/Grid.h"
#include "FEM/FEM.h"
#include <iostream>
#include <fstream>

int main()
{
	Grid grid;

	grid.ReadGridJSON();
	grid.CreateGrid();

	auto numEdge = grid.GetNumberEdge(Grid::Edge{ {275, 276, 282}, {0, -1} });
	
	std::cout << numEdge;
	std::cout << std::endl << std::endl;

	auto numsEdges = grid.GetNumbersEdges(915);

	for(const auto& num: numsEdges)
		std::cout << num << " ";
	std::cout << std::endl << std::endl;

	Grid::Edge edge = grid.GetEdge(1777);

	for (const auto& v : edge.vertexes)
		std::cout << v << " ";
	std::cout << std::endl;
	for (const auto& n : edge.elemNums)
		std::cout << n << " ";
	std::cout << std::endl << std::endl;

	edge = grid.GetEdge(2044);

	for (const auto& v : edge.vertexes)
		std::cout << v << " ";
	std::cout << std::endl;
	for (const auto& n : edge.elemNums)
		std::cout << n << " ";
	std::cout << std::endl << std::endl;

	edge = grid.GetEdge(2045);

	for (const auto& v : edge.vertexes)
		std::cout << v << " ";
	std::cout << std::endl;
	for (const auto& n : edge.elemNums)
		std::cout << n << " ";
	std::cout << std::endl << std::endl;

	edge = grid.GetEdge(2046);

	for (const auto& v : edge.vertexes)
		std::cout << v << " ";
	std::cout << std::endl;
	for (const auto& n : edge.elemNums)
		std::cout << n << " ";
	std::cout << std::endl << std::endl;


	grid.WriteGrid();

	FEM fem(grid, [](const std::array<double, SIZE_NODE>& args)
		{auto&& [x, y, z] = args; return x + y + z; });
	
	fem.ReadParametersJSON();
	fem.CollectSLAE();
	fem.ConsiderBoundaryConditions();
	fem.Solve();
	//fem.CheckSolution(std::cout);
	fem.WriteSolution();

	return 0;
}

