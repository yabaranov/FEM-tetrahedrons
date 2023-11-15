#include "LocalAssembly.h"
#include "FEM.h"
#include "differentiation.h"

void LocalAssembly::GetFiniteVectorBC_2(FiniteVector<SIZE_EDGE>& lV, int edgeNum)
{
	double areaBoundaryEdge = GetAreaBoundaryEdge(Grid::TYPE_BOUNDARY_CONDITION::SECOND, edgeNum);
	FiniteVector<SIZE_EDGE> tetta = GetInterpolantTetta(edgeNum);	
	lV = areaBoundaryEdge / 12 * m_A * tetta;
}

void LocalAssembly::GetFiniteMatrixVectorBC_3(FiniteMatrix<SIZE_EDGE>& lM, FiniteVector<SIZE_EDGE>& lV, int edgeNum)
{
	double areaBoundaryEdge = GetAreaBoundaryEdge(Grid::TYPE_BOUNDARY_CONDITION::THIRD, edgeNum);
	lM = m_fem.Betta(m_grid.GetElement(m_grid.GetBoundaryEdgeBC_3(edgeNum).elemNum).numberFormula) * areaBoundaryEdge / 12 * m_A;
	FiniteVector<SIZE_EDGE> u_betta = GetInterpolantU_betta(edgeNum);	
	lV = lM * u_betta;
}

void LocalAssembly::GetFiniteMatrixVectorElement(FiniteMatrix<SIZE_ELEMENT>& lM, FiniteVector<SIZE_ELEMENT>& lV, int elemNum)
{
	FiniteMatrix<SIZE_ELEMENT> G;
	double detD = GetG(G, elemNum);

	lM = m_fem.Lambda(m_grid.GetElement(elemNum).numberFormula) * detD / 6 * G +
		m_fem.Gamma(m_grid.GetElement(elemNum).numberFormula) * detD / 120 * m_C;
	
	FiniteVector<SIZE_ELEMENT> f = GetInterpolantF(elemNum);	
	lV = detD / 120 * m_C * f;
}

double LocalAssembly::GetG(FiniteMatrix<SIZE_ELEMENT>& G, int elemNum)
{
	FiniteMatrix<SIZE_ELEMENT> D = GetD(elemNum);
	
	double detD = std::abs(D.InverseMatrix());

	for (int k = 0; k < G.size(); k++)
		for (int j = 0; j < G.size(); j++)
			G[k][j] = D[k][1] * D[j][1] + D[k][2] * D[j][2] + D[k][3] * D[j][3];

	return detD;
}

FiniteMatrix<SIZE_ELEMENT> LocalAssembly::GetD(int elemNum)
{
	FiniteMatrix<SIZE_ELEMENT> D;
	for (int k = 0; k < D.size(); k++)
		D[0][k] = 1;
	for (int j = 0; j < D.size(); j++)
		for (int k = 1; k < D.size(); k++)
			D[k][j] = m_grid.GetNode(m_grid.GetElement(elemNum).vertexes[j])[k - 1];

	return D;
}

FiniteVector<SIZE_EDGE> LocalAssembly::GetInterpolantTetta(int edgeNum)
{
	FiniteVector<SIZE_NODE> n = Getn(Grid::TYPE_BOUNDARY_CONDITION::SECOND, edgeNum);
	FiniteVector<SIZE_EDGE> tetta;
	for (int j = 0; j < tetta.size(); j++)	
		tetta[j] = m_fem.NumericalTetta(n, edgeNum, j);
	return tetta;
}

FiniteVector<SIZE_EDGE> LocalAssembly::GetInterpolantU_betta(int edgeNum)
{
	FiniteVector<SIZE_NODE> n = Getn(Grid::TYPE_BOUNDARY_CONDITION::THIRD, edgeNum);
	FiniteVector<SIZE_EDGE> u_betta;
	for (int j = 0; j < u_betta.size(); j++)
		u_betta[j] = m_fem.NumericalU_betta(n, edgeNum, j);
	return u_betta;
}

FiniteVector<SIZE_ELEMENT> LocalAssembly::GetInterpolantF(int elemNum)
{
	FiniteVector<SIZE_ELEMENT> f;
	for (int j = 0; j < f.size(); j++)
		f[j] = m_fem.NumericalF(elemNum, j);
	return f;
}

double LocalAssembly::GetAreaBoundaryEdge(Grid::TYPE_BOUNDARY_CONDITION typeBoundaryCondition, int edgeNum)
{
	const Grid::BoundaryEdge& BC = (typeBoundaryCondition == Grid::TYPE_BOUNDARY_CONDITION::SECOND) ? m_grid.GetBoundaryEdgeBC_2(edgeNum) : m_grid.GetBoundaryEdgeBC_3(edgeNum);

	FiniteMatrix<SIZE_NODE> coords;

	for (int j = 0; j < coords.size(); j++)
		for (int k = 0; k < coords.size(); k++)
			coords[j][k] = m_grid.GetNode(BC.vertexes[j])[k];
			
	double px = (coords[1][1] - coords[0][1]) * (coords[2][2] - coords[0][2]) - (coords[2][1] - coords[0][1]) * (coords[1][2] - coords[0][2]);
	double py = (coords[2][0] - coords[0][0]) * (coords[1][2] - coords[0][2]) - (coords[1][0] - coords[0][0]) * (coords[2][2] - coords[0][2]);
	double pz = (coords[1][0] - coords[0][0]) * (coords[2][1] - coords[0][1]) - (coords[2][0] - coords[0][0]) * (coords[1][1] - coords[0][1]);

	return std::sqrt(px * px + py * py + pz * pz) / 2;
}

FiniteVector<SIZE_NODE> LocalAssembly::Getn(Grid::TYPE_BOUNDARY_CONDITION typeBoundaryCondition, int edgeNum)
{
	const Grid::BoundaryEdge& BC = (typeBoundaryCondition == Grid::TYPE_BOUNDARY_CONDITION::SECOND) ? m_grid.GetBoundaryEdgeBC_2(edgeNum) : m_grid.GetBoundaryEdgeBC_3(edgeNum);

	FiniteVector<SIZE_NODE> n = vectorProduct
	(
		FiniteVector<SIZE_NODE>({ m_grid.GetNode(BC.vertexes[1]) }) - FiniteVector<SIZE_NODE>({ m_grid.GetNode(BC.vertexes[0]) }),
		FiniteVector<SIZE_NODE>({ m_grid.GetNode(BC.vertexes[2]) }) - FiniteVector<SIZE_NODE>({ m_grid.GetNode(BC.vertexes[0]) })
	);

	n /= sqrt(n * n);

	int k;
	for (k = 0; k < SIZE_ELEMENT; k++)
		if (binarySearch(BC.vertexes, m_grid.GetElement(BC.elemNum).vertexes[k], 0, BC.vertexes.size() - 1) == -1)
			break;

	FiniteVector<SIZE_NODE> vectorOutBoundaryEdge = FiniteVector<SIZE_NODE>({ m_grid.GetNode(m_grid.GetElement(BC.elemNum).vertexes[k]) }) - FiniteVector<SIZE_NODE>({ m_grid.GetNode(BC.vertexes[0]) });

	if (n * vectorOutBoundaryEdge > 0)
		n = -n;

	return n;
}