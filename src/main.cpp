#include "grid/Grid.h"
#include "FEM/FEM.h"
#include <iostream>
#include <fstream>

int main()
{
	Grid grid;

	grid.ReadGridJSON();
	grid.CreateGrid();
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

