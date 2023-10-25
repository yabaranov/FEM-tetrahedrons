#include "SparseSLAE.h"
#include "../functions.h"

void SparseSLAE::RegenerationProfileFormat()
{
    std::vector<int> igp(m_M.di.size() + 1);

    igp[0] = 0;
    for (int i = 0; i < m_M.di.size(); i++)
        (m_M.ig[i + 1] > m_M.ig[i]) ? igp[i + 1] = igp[i] + i - m_M.jg[m_M.ig[i]] : igp[i + 1] = igp[i];

    for (int i = 0; i < m_M.di.size(); i++)
    {
        int column = i - (igp[i + 1] - igp[i]);
        int s_point = m_M.ig[i];
        for (int j = igp[i]; j < igp[i + 1]; j++, column++)
            if (column != m_M.jg[s_point])
            {
                m_M.ggl.insert(m_M.ggl.begin() + j, 1, 0);
                m_M.ggu.insert(m_M.ggu.begin() + j, 1, 0);
            }                    
            else s_point++;
    }

    m_M.ig = igp;
 
}

void SparseSLAE::LUDecomposition()
{
    for (int i = 0; i < m_M.di.size(); i++)
    {
        int j = i - (m_M.ig[i + 1] - m_M.ig[i]);
        double sd = 0;
        for (int m = m_M.ig[i]; m < m_M.ig[i + 1]; m++, j++)
        {           
            int mi = m_M.ig[i];
            int mj = m_M.ig[j];

            int kol_r = m - m_M.ig[i] - m_M.ig[j + 1] + m_M.ig[j];
            (kol_r < 0) ? mj -= kol_r : mi += kol_r;

            double sl = 0;
            double su = 0;

            for (mi = mi; mi < m; mi++, mj++)
            {
                sl += m_M.ggl[mi] * m_M.ggu[mj];
                su += m_M.ggu[mi] * m_M.ggl[mj];
            }

            m_M.ggl[m] -= sl;
            m_M.ggu[m] = (m_M.ggu[m] - su) / m_M.di[j];
            sd += m_M.ggl[m] * m_M.ggu[m];
        }

        m_M.di[i] -= sd;
    }
}

void SparseSLAE::LUSolve()
{  
    RegenerationProfileFormat();

    LUDecomposition();

    for (int i = 0; i < m_M.di.size(); i++)
    {
        double sum = 0;
        int vect_iter = i - (m_M.ig[i + 1] - m_M.ig[i]);
        for (int j = m_M.ig[i]; j < m_M.ig[i + 1]; j++, vect_iter++)
            sum += m_M.ggl[j] * m_x[vect_iter];

        m_x[i] = (m_b[i] - sum) / m_M.di[i];
    }

    for (int i = m_M.di.size() - 1; i >= 0; i--)
    {
        double sum = 0;
        int vect_iter = i - (m_M.ig[i + 1] - m_M.ig[i]);
        for (int j = m_M.ig[i]; j < m_M.ig[i + 1]; j++, vect_iter++)
            m_x[vect_iter] -= m_M.ggu[j] * m_x[i];;
    }

}

int SparseSLAE::BÑGSTABLU(double eps, int maxIterations)
{
    int n = m_M.di.size();

    Vector r(n);
    Vector r_0(n);
    Vector p(n);
    Vector s(n);
    Vector LAUp(n);
    Vector LAUs(n);

    Vector L(m_M.ig[n]);
    Vector U(m_M.ig[n]);
    Vector D(n);

    double alpha, omega, betta, dotPrr_0r;
    double normb = sqrt(m_b * m_b);
    double residual;
    IncompleteLUDecomposition(L, U, D);

    StraightMove(L, D, m_b, r_0);
    r = r_0;
    BackMove(U, r_0, p);

    int k = 1;
    for (; k < maxIterations && (residual = sqrt(r * r) / normb) > eps; k++)
    {
 
        BackMove(U, p, LAUp);
        LAUp = m_M * LAUp;
        StraightMove(L, D, LAUp, LAUp);
        dotPrr_0r = r_0 * r;
        alpha = dotPrr_0r / (LAUp * r_0);

        s = r - alpha * LAUp;

        BackMove(U, s, LAUs);
        LAUs = m_M * LAUs;
        StraightMove(L, D, LAUs, LAUs);
        omega = (LAUs * s) / (LAUs * LAUs);

        m_x += alpha * p + omega * s;

        r = s - omega * LAUs;

        betta = alpha * (r * r_0) / (omega * dotPrr_0r);

        p = r + betta * (p - omega * LAUp);

    }

    BackMove(U, m_x, m_x);

    return k;
}

int SparseSLAE::LosLU(double eps, int maxIterations)
{
    int n = m_M.di.size();

    Vector r(n);
    Vector z(n);
    Vector p(n);
    Vector LAUr(n);
    Vector Ur(n);

    Vector L(m_M.ig[n]);
    Vector U(m_M.ig[n]);
    Vector D(n);

    double alpha, betta, dotPrpp;
    double normb = sqrt(m_b * m_b);
    double residual;
    IncompleteLUDecomposition(L, U, D);

    StraightMove(L, D, m_b, r);
    BackMove(U, r, z);
    StraightMove(L, D, m_M * z, p);

    int k = 1;
    for (; k < maxIterations && (residual = sqrt(r * r) / normb) > eps; k++)
    {
        dotPrpp = p * p;

        alpha = (p * r) / dotPrpp;

        m_x += alpha * z;
        r -= alpha * p;

        BackMove(U, r, Ur);
        LAUr = m_M * Ur;
        StraightMove(L, D, LAUr, LAUr);

        betta = (-p * LAUr) / dotPrpp;

        z = Ur + betta * z;
        p = LAUr + betta * p;

    }

    return k;
}

void SparseSLAE::StraightMove(const Vector& L, const Vector& D, const Vector& v1, Vector& v2) const
{
    int n = m_M.di.size();
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = m_M.ig[i]; j < m_M.ig[i + 1]; j++)
            sum += v2[m_M.jg[j]] * L[j];

        v2[i] = (v1[i] - sum) / D[i];
    }
}

void SparseSLAE::BackMove(const Vector& U, const Vector& v1, Vector& v2) const
{
    int n = m_M.di.size();
    for (int i = 0; i < n; i++)
        v2[i] = v1[i];

    for (int i = n - 1; i >= 0; i--)       
        for (int j = m_M.ig[i]; j < m_M.ig[i + 1]; j++)
            v2[m_M.jg[j]] -= v2[i] * U[j];
       
}

void SparseSLAE::IncompleteLUDecomposition(Vector& L, Vector& U, Vector& D) const
{
    int n = m_M.di.size();

    D = m_M.di;
    L = m_M.ggl;
    U = m_M.ggu;

    for (int i = 0; i < n; i++)
    {
        double d = 0;
        int temp = m_M.ig[i];
        for (int j = m_M.ig[i]; j < m_M.ig[i + 1]; j++)
        {
            double ls = 0;
            double us = 0;
            for (int h = m_M.ig[m_M.jg[j]], k = temp; h < m_M.ig[m_M.jg[j] + 1] && k < j;)
                if (m_M.jg[k] == m_M.jg[h])
                {
                    ls += L[k] * U[h];
                    us += L[h++] * U[k++];
                }
                else (m_M.jg[k] < m_M.jg[h]) ? k++ : h++;

            L[j] -= ls;
            U[j] = (U[j] - us) / D[m_M.jg[j]];
            d += L[j] * U[j];
        }
        D[i] -= d;
    }
}

Vector SparseSLAE::SparseMatrix::operator*(const Vector& v) const
{
	Vector result(v.size());

	for (int i = 0; i < v.size(); i++)
	{
		result[i] = di[i] * v[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			result[i] += ggl[k] * v[j];
			result[j] += ggu[k] * v[i];
		}
	}		

    return result;
}

void SparseSLAE::ZeroingRow(int row)
{
    int ibeg = m_M.ig[row];
    int iend = m_M.ig[row + 1];
    for (int i = ibeg; i < iend; i++)
        m_M.ggl[i] = 0;

    for (int j = row + 1; j < m_b.size(); j++)
    {
        int jbeg = m_M.ig[j];
        int jend = m_M.ig[j + 1];
        int index = binarySearch(m_M.jg, row, jbeg, jend - 1);
        if (index != -1)
            m_M.ggu[index] = 0;  
    }
}

void SparseSLAE::ZeroingSLAE()
{
    m_M.ggl = 0.0;
    m_M.ggu = 0.0;
    m_M.di = 0.0;
    m_b = 0.0;
    m_x = 0.0;   
}

void SparseSLAE::GenerateSLAE(const std::vector<Grid::Element>& elements, int N, int neighboringNodes)
{
    int memory = N * neighboringNodes;
    std::vector<std::vector<int>> list(2, std::vector<int>(memory, 0));

    std::vector<int>& ig = m_M.ig;
    ig.reserve(N);
    std::vector<int>& jg = m_M.jg;
    jg.reserve(memory);
    std::vector<int> listbeg(N, 0);
    int listSize = 0;

    for (const auto& elem : elements)
    {
        for (int i = 0; i < elem.vertexes.size(); i++)
        {          
            int k = elem.vertexes[i];

            for (int j = i + 1; j < elem.vertexes.size(); j++)
            {
                int ind1 = k;

                int ind2 = elem.vertexes[j];
                if (ind2 < ind1)
                {
                    ind1 = ind2;
                    ind2 = k;
                }
                int iaddr = listbeg[ind2];
                if (!iaddr)
                {
                    listSize++;
                    listbeg[ind2] = listSize;
                    list[0][listSize] = ind1;
                    list[1][listSize] = 0;
                }
                else
                {
                    while (list[0][iaddr] < ind1 && list[1][iaddr] > 0)
                    {
                        iaddr = list[1][iaddr];
                    }
                    if (list[0][iaddr] > ind1)
                    {
                        listSize++;
                        list[0][listSize] = list[0][iaddr];
                        list[1][listSize] = list[1][iaddr];
                        list[0][iaddr] = ind1;
                        list[1][iaddr] = listSize;
                    }
                    else
                    {
                        if (list[0][iaddr] < ind1)
                        {
                            listSize++;
                            list[1][iaddr] = listSize;
                            list[0][listSize] = ind1;
                            list[1][listSize] = 0;
                        }
                    }
                }
            }
        }
    }

    ig.push_back(0);
    for (int i = 0; i < N; i++)
    {
        ig.push_back(ig[i]);
        int iaddr = listbeg[i];
        while (iaddr != 0)
        {
            jg.push_back(list[0][iaddr]);
            ig[i + 1]++;
            iaddr = list[1][iaddr];
        }
    }

    m_M.ggl.resize(m_M.jg.size());
    m_M.ggu.resize(m_M.jg.size());
    m_M.di.resize(N);
    m_b.resize(N); 
    m_x.resize(N);
}

int SparseSLAE::Solve(TYPE_SOLVER typeSolver)
{
    switch (typeSolver)
    {
    case TYPE_SOLVER::LU:
        LUSolve(); return 0;
        break;
    case TYPE_SOLVER::LOSLU:
        return LosLU(1e-15, 10000);
        break;
    case TYPE_SOLVER::BÑGSTABLU:
        return BÑGSTABLU(1e-15, 10000);
        break;
    }
}
