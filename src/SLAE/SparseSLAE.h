#ifndef SPARSE_SLAE
#define SPARSE_SLAE
#include "../vector/MatrixVector.h"
#include "../grid/Grid.h"

class SparseSLAE
{
public:	

	struct SparseMatrix
	{
		Vector di;
		std::vector<int> ig;
		std::vector<int> jg;
		Vector ggl;
		Vector ggu;
		Vector operator*(const Vector& v) const;
	};

	enum class TYPE_SOLVER { LU, LOSLU, BÑGSTABLU };
	SparseSLAE() {}
	
	void GenerateSLAE(const std::vector<Grid::Element>& elements, int N, int neighbors);

	int Solve(TYPE_SOLVER typeSolver);
	void ZeroingRow(int row);
	void ZeroingSLAE();

	SparseMatrix& Set_M() { return m_M; }
	const SparseMatrix& Get_M() const { return m_M; }
	Vector& Set_b() { return m_b; }
	const Vector& Get_b() const { return m_b; }
	Vector& Set_x() { return m_x; }
	const Vector& Get_x() const { return m_x; }


	~SparseSLAE() {}

private:

	void LUSolve();
	int LosLU(double eps, int maxIterations);
	int BÑGSTABLU(double eps, int maxIterations);

	SparseMatrix m_M;
	Vector m_b;
	Vector m_x;

	void RegenerationProfileFormat();
	void BackMove(const Vector& U, const Vector& v1, Vector& v2) const;
	void StraightMove(const Vector& L, const Vector& D, const Vector& v1, Vector& v2) const;
	void LUDecomposition();
	void IncompleteLUDecomposition(Vector& L, Vector& U, Vector& D) const;
};
#endif
