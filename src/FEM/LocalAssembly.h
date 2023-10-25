#ifndef LOCAL_ASSEMBLY
#define LOCAL_ASSEMBLY
#include "FEM.h"
#include "../grid/Grid.h"
#include "../vector/FiniteMatrixVector.h"
#include "../functions.h"

class LocalAssembly
{
public:

	LocalAssembly(const FEM& fem) : m_fem(fem), m_grid(fem.GetGrid()) {}
	void GetFiniteMatrixVectorElement(FiniteMatrix<SIZE_ELEMENT>& lM, FiniteVector<SIZE_ELEMENT>& lV, int elemNum);
	void GetFiniteVectorBC_2(FiniteVector<SIZE_EDGE>& lV, int edgeNum);
	void GetFiniteMatrixVectorBC_3(FiniteMatrix<SIZE_EDGE>& lM, FiniteVector<SIZE_EDGE>& lV, int edgeNum);

private:
	enum TYPE_BOUNDARY_CONDITION {SECOND, THIRD};

	const FEM& m_fem;
	const Grid& m_grid;

	const FiniteMatrix<SIZE_ELEMENT> m_C =
	{{
		FiniteVector<SIZE_ELEMENT>({ 2, 1, 1, 1 }),		
		FiniteVector<SIZE_ELEMENT>({ 1, 2, 1, 1 }),			
		FiniteVector<SIZE_ELEMENT>({ 1, 1, 2, 1 }),			
		FiniteVector<SIZE_ELEMENT>({ 1, 1, 1, 2 }),			
	}};

	const FiniteMatrix<SIZE_EDGE> m_A =
	{{
		FiniteVector<SIZE_EDGE>({ 2, 1, 1 }),		
		FiniteVector<SIZE_EDGE>({ 1, 2, 1 }),		
		FiniteVector<SIZE_EDGE>({ 1, 1, 2 }),		
	}};

	double GetG(FiniteMatrix<SIZE_ELEMENT>& G, int elemNum);
	FiniteMatrix<SIZE_ELEMENT> GetD(int elemNum);
	double GetAreaEdge(TYPE_BOUNDARY_CONDITION typeBoundaryCondition, int edgeNum);

	FiniteVector<SIZE_NODE> Getn(TYPE_BOUNDARY_CONDITION typeBoundaryCondition, int edgeNum);

	FiniteVector<SIZE_ELEMENT> GetInterpolantF(int elemNum);
	FiniteVector<SIZE_EDGE> GetInterpolantTetta(int edgeNum);
	FiniteVector<SIZE_EDGE> GetInterpolantU_betta(int edgeNum);

};

#endif