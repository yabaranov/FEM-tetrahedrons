#include <iomanip>
#include <rapidjson/document.h>
#include <rapidjson/error/en.h>
#include <fstream>
#include "FEM.h"
#include "LocalAssembly.h"
#include "differentiation.h"

bool FEM::ReadParametersJSON()
{
	std::string JSONPath = "../res/input/parameters.json";
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

	auto lambdasIt = document.FindMember("lambdas");
	if (lambdasIt != document.MemberEnd())
	{
		const auto lambdasArray = lambdasIt->value.GetArray();
		m_lambdas.resize(lambdasArray.Size());
		for (size_t i = 0; i < lambdasArray.Size(); i++)
			m_lambdas[i] = lambdasArray[i].GetDouble();
	}

	auto gammasIt = document.FindMember("gammas");
	if (gammasIt != document.MemberEnd())
	{
		const auto gammasArray = gammasIt->value.GetArray();
		m_gammas.resize(gammasArray.Size());
		for (size_t i = 0; i < gammasArray.Size(); i++)
			m_gammas[i] = gammasArray[i].GetDouble();
	}

	auto bettasIt = document.FindMember("bettas");
	if (bettasIt != document.MemberEnd())
	{
		const auto bettasArray = bettasIt->value.GetArray();
		m_bettas.resize(bettasArray.Size());
		for (size_t i = 0; i < bettasArray.Size(); i++)
			m_bettas[i] = bettasArray[i].GetDouble();
	}

	auto typeSolverIt = document.FindMember("typeSolver");
	if (typeSolverIt != document.MemberEnd())
	{
		std::string typeSolver = typeSolverIt->value.GetString();

		if (typeSolver == "BCGSTABLU")
			m_typeSolver = SparseSLAE::TYPE_SOLVER::BCGSTABLU;
		else if(typeSolver == "LOSLU")
			m_typeSolver = SparseSLAE::TYPE_SOLVER::LOSLU;
		else if(typeSolver == "LU")
			m_typeSolver = SparseSLAE::TYPE_SOLVER::LU;
		else
		{
			std::cerr << "Error type solver!"; 
			return false;
		}
	}

	return true;
}

void FEM::CollectSLAE()
{
	m_sparseSLAE.GenerateSLAE(m_grid.Elements(), m_grid.SizeNodes(), m_grid.GetNumberNeighboringNodes());

	LocalAssembly lA(*this);
	for (int elemNum = 0; elemNum < m_grid.SizeElements(); elemNum++)
	{
		FiniteMatrix<SIZE_ELEMENT> lM;
		FiniteVector<SIZE_ELEMENT> lV;
		lA.GetFiniteMatrixVectorElement(lM, lV, elemNum);

		AddFiniteMatrix(lM, m_grid.GetElement(elemNum).vertexes);
		AddFiniteVector(lV, m_grid.GetElement(elemNum).vertexes);
	}
}

void FEM::ConsiderBoundaryConditions()
{
	ConsiderBC_2();
	ConsiderBC_3();
	ConsiderBC_1();
}

void FEM::ConsiderBC_2()
{
	LocalAssembly lA(*this);
	
	for (int edgeNum = 0; edgeNum < m_grid.SizeBoundaryEdgesBC_2(); edgeNum++)
	{
		FiniteVector<SIZE_EDGE> lV;
		lA.GetFiniteVectorBC_2(lV, edgeNum);
		AddFiniteVector(lV, m_grid.GetBoundaryEdgeBC_2(edgeNum).vertexes);
	}
}

void FEM::ConsiderBC_3()
{
	LocalAssembly lA(*this);
	
	for (int edgeNum = 0; edgeNum < m_grid.SizeBoundaryEdgesBC_3(); edgeNum++)
	{
		FiniteMatrix<SIZE_EDGE> lM;
		FiniteVector<SIZE_EDGE> lV;
		lA.GetFiniteMatrixVectorBC_3(lM, lV, edgeNum);

		AddFiniteMatrix(lM, m_grid.GetBoundaryEdgeBC_3(edgeNum).vertexes);
		AddFiniteVector(lV, m_grid.GetBoundaryEdgeBC_3(edgeNum).vertexes);
	}
}

void FEM::ConsiderBC_1()
{
	SparseSLAE::SparseMatrix& M = m_sparseSLAE.Set_M();
	Vector& b = m_sparseSLAE.Set_b();

	for (int nodeNum = 0; nodeNum < m_grid.SizeNodesBC_1(); nodeNum++)
	{
		int numberNodeBC_1 = m_grid.GetNumberNodeBC_1(nodeNum);

		M.di[numberNodeBC_1] = 1;
		m_sparseSLAE.ZeroingRow(numberNodeBC_1);
		b[numberNodeBC_1] = m_u_g(m_grid.GetNode(numberNodeBC_1));
	}
}

double FEM::NumericalF(int elemNum, int verNum) const
{
	return  -Lambda(m_grid.GetElement(elemNum).numberFormula) * operatorLaplace(m_u_g, m_grid.GetNode(m_grid.GetElement(elemNum).vertexes[verNum]))
		+ Gamma(m_grid.GetElement(elemNum).numberFormula) * m_u_g(m_grid.GetNode(m_grid.GetElement(elemNum).vertexes[verNum]));
}

double FEM::NumericalTetta(const FiniteVector<SIZE_NODE>& n, int edgeNum, int verNum) const
{
	return Lambda(m_grid.GetElement(m_grid.GetBoundaryEdgeBC_2(edgeNum).elemNum).numberFormula) * gradient(m_u_g, m_grid.GetNode(m_grid.GetBoundaryEdgeBC_2(edgeNum).vertexes[verNum])) * n;
}

double FEM::NumericalU_betta(const FiniteVector<SIZE_NODE>& n, int edgeNum, int verNum) const
{
	return Lambda(m_grid.GetElement(m_grid.GetBoundaryEdgeBC_3(edgeNum).elemNum).numberFormula) / Betta(m_grid.GetElement(m_grid.GetBoundaryEdgeBC_3(edgeNum).elemNum).numberFormula) *
		gradient(m_u_g, m_grid.GetNode(m_grid.GetBoundaryEdgeBC_3(edgeNum).vertexes[verNum])) * n + m_u_g(m_grid.GetNode(m_grid.GetBoundaryEdgeBC_3(edgeNum).vertexes[verNum]));
}

std::ostream& FEM::CheckSolution(std::ostream& os) const
{
	const Vector& x = m_sparseSLAE.Get_x();

	double sum_1 = 0, sum_2 = 0;
	
	os << "+-------+-------+-------+----------------------+----------------------+----------------------+\n";
	os << "|   X   |   Y   |   Z   |          u           |          u*          |         |u-u*|       |\n";
	os << "+-------+-------+-------+----------------------+----------------------+----------------------+\n";
	
	
	for (int i = 0; i < m_grid.SizeNodes(); i++)
	{
		auto node = m_grid.GetNode(i);
		for (int j = 0; j < node.size(); j++)
			os << "|" << std::fixed << std::setw(7) << std::setprecision(3) << node[j];
		
		os << "|" << std::scientific << std::setw(22) << std::setprecision(DBL_DIG) << x[i] << "|"
			<< std::scientific << std::setw(22) << std::setprecision(DBL_DIG) << m_u_g(node) << "|"
			<< std::scientific << std::setw(22) << std::setprecision(DBL_DIG) << std::abs(x[i] - m_u_g(node)) << "|\n";

		if (i % m_grid.GetGridNesting() == 0)
		{
			sum_1 += (x[i] - m_u_g(node)) * (x[i] - m_u_g(node));
			sum_2 += m_u_g(node) * m_u_g(node);
		}
		
	}
	
	os << "+--------------------------------------------------------------------------------------------+\n";
	os << "|   ||u-u*||/||u*||    " << std::scientific << std::setw(70) << std::setprecision(DBL_DIG) << std::sqrt(sum_1 / sum_2) << "|\n";
	os << "+--------------------------------------------------------------------------------------------+\n";
	
	return os;
}

bool FEM::WriteSolution() const
{
	std::string filePath = "../res/output/solution.txt";
	std::ofstream solution;
	solution.open(filePath.c_str(), std::ios::out | std::ios::binary);
	if (!solution.is_open())
	{
		std::cerr << "Failed to open file: " << filePath << std::endl;
		return false;
	}

	const Vector& x = m_sparseSLAE.Get_x();

	for (int i = 0; i < m_grid.SizeNodes(); i++)
	{
		auto node = m_grid.GetNode(i);
		for (int j = 0; j < node.size(); j++)
			solution << " " << std::fixed << std::setw(7) << std::setprecision(3) << node[j];

		solution << " " << std::scientific << std::setw(22) << std::setprecision(DBL_DIG) << x[i] << "\n";
	}

	return true;
}

