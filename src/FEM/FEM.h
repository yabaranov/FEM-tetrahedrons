#pragma once
#include "../SLAE/SparseSLAE.h"
#include "../grid/Grid.h"
#include "../vector/FiniteMatrixVector.h"
#include "../functions.h"


class FEM
{
public:

	FEM(const Grid& grid, std::function<double(const std::array<double, SIZE_NODE>& args)> u_g) : m_grid(grid), m_u_g(u_g){}
	~FEM() {}
	bool ReadParametersJSON();
	void CollectSLAE();
	void Solve() { m_sparseSLAE.Solve(m_typeSolver); }
	void ConsiderBoundaryConditions();
	std::ostream& Output(std::ostream& os) const;

	const Grid& GetGrid() const { return m_grid; }

	double NumericalF(int elemNum, int verNum) const;
	double NumericalTetta(FiniteVector<SIZE_NODE> n, int edgeNum, int verNum) const;
	double NumericalU_betta(FiniteVector<SIZE_NODE> n, int edgeNum, int verNum) const;
	double Lambda(int numberFormula) const { return m_lambdas[numberFormula]; }
	double Betta(int numberFormula) const { return m_bettas[numberFormula]; }
	double Gamma(int numberFormula) const { return m_gammas[numberFormula]; }

private:

	SparseSLAE::TYPE_SOLVER m_typeSolver = SparseSLAE::TYPE_SOLVER::BÑGSTABLU;
	SparseSLAE m_sparseSLAE;
	const Grid& m_grid;
	
	template<int n>
	void AddFiniteMatrix(const FiniteMatrix<n>& lM, const std::array<int, n>& L);

	template<int n>
	void AddFiniteVector(const FiniteVector<n>& lV, const std::array<int, n>& L);

	void ConsiderBC_1();
	void ConsiderBC_2();
	void ConsiderBC_3();

	std::vector<double> m_lambdas;	
	std::vector<double> m_gammas;	
	std::vector<double> m_bettas;

	std::function<double(const std::array<double, SIZE_NODE>& args)> m_u_g;	
};

template<int n>
void FEM::AddFiniteVector(const FiniteVector<n>& lV, const std::array<int, n>& L)
{
	Vector& b = m_sparseSLAE.Set_b();
	for (int i = 0; i < n; i++)
		b[L[i]] += lV[i];
}

template<int n>
void FEM::AddFiniteMatrix(const FiniteMatrix<n>& lM, const std::array<int, n>& L)
{
	SparseSLAE::SparseMatrix& M = m_sparseSLAE.Set_M();

	for (int i = 0; i < n; i++)
	{
		M.di[L[i]] += lM[i][i];

		int ibeg = M.ig[L[i]];
		int iend = M.ig[L[i] + 1];
		for (int j = 0; j < i; j++)
		{		
			int index = binarySearch(M.jg, L[j], ibeg, iend - 1);

			M.ggu[index] += lM[j][i];
			M.ggl[index] += lM[i][j];
			ibeg++;
		}
	}
}