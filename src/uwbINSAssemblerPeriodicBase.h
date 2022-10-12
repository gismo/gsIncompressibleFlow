/** @file uwbINSAssemblerPeriodicBase.h

Author(s): H. Hornikova
*/

#pragma once

namespace gismo
{

template<class T>
class uwbINSAssemblerPeriodicBase
{

public:
    uwbINSAssemblerPeriodicBase(uwbINSSolverParams<T>& params, uwbINSBlockAssembler<T>& blockAssembler)
        : m_blockAssemblerRef(blockAssembler)
    {
        m_transformMatrix = params.getBCs().getTransformMatrix();

        m_invTransformMatrix = m_transformMatrix.transpose();

        initMembers();
    }

    virtual ~uwbINSAssemblerPeriodicBase()
    {
    }

protected:

    void initMembers()
    {
        initPeriodic();

        m_perSolution.setZero(m_dofs, 1);

        m_bNonZerosExact = false;
    }

    virtual void reinitMembers() { initMembers(); }

    //initialize periodic stuff and mappings
    void initPeriodic()
    {
        const gsDofMapper & Umapper = m_blockAssemblerRef.getMappers().front();
        const gsDofMapper & Pmapper = m_blockAssemblerRef.getMappers().back();
        int blockUdofs = m_blockAssemblerRef.getUdofs();
        int blockPdofs = m_blockAssemblerRef.getPdofs();

        m_uIsEliminated.setZero(blockUdofs, 1);
        m_pIsEliminated.setZero(blockPdofs, 1);
        m_perDofs1_u.clear();
        m_perDofs2_u.clear();
        m_perDofs1_p.clear();
        m_perDofs2_p.clear();

        for (typename gsBoundaryConditions<T>::const_ppiterator it = m_blockAssemblerRef.getBCs().periodicBegin(); it != m_blockAssemblerRef.getBCs().periodicEnd(); ++it)
        {
            index_t patch1 = it->first().patch;
            index_t patch2 = it->second().patch;

            // local indices of periodic dofs on both sides
            gsMatrix<index_t> pd1_u, pd2_u, pd1_p, pd2_p;
            m_blockAssemblerRef.getBases().front().basis(patch1).matchWith(*it, m_blockAssemblerRef.getBases().front().basis(patch2), pd1_u, pd2_u);
            m_blockAssemblerRef.getBases().back().basis(patch1).matchWith(*it, m_blockAssemblerRef.getBases().back().basis(patch2), pd1_p, pd2_p);

            // global indices of the periodic dofs
            gsMatrix<index_t> gl_pd1_u, gl_pd2_u, gl_pd1_p, gl_pd2_p;
            Umapper.localToGlobal(pd1_u, patch1, gl_pd1_u);
            Umapper.localToGlobal(pd2_u, patch2, gl_pd2_u);
            Pmapper.localToGlobal(pd1_p, patch1, gl_pd1_p);
            Pmapper.localToGlobal(pd2_p, patch2, gl_pd2_p);

            // append the global indices to vectors of periodic dofs
            for (index_t i = 0; i < gl_pd1_u.rows(); i++)
            {
                if ((Umapper.is_free_index(gl_pd1_u(i))) && (!m_uIsEliminated(gl_pd2_u(i))) && (!m_uIsEliminated(gl_pd1_u(i)))) {
                    m_perDofs1_u.push_back(gl_pd1_u(i));
                    m_perDofs2_u.push_back(gl_pd2_u(i));
                    m_uIsEliminated(gl_pd2_u(i)) = true;
                }
            }

            for (index_t i = 0; i < gl_pd1_p.rows(); i++)
            {
                if ((Pmapper.is_free_index(gl_pd1_p(i)) && (!m_pIsEliminated(gl_pd2_p(i)))) && (!m_pIsEliminated(gl_pd1_p(i)))) {
                    m_perDofs1_p.push_back(gl_pd1_p(i));
                    m_perDofs2_p.push_back(gl_pd2_p(i));
                    m_pIsEliminated(gl_pd2_p(i)) = true;
                }
            }
        }

        m_perUdofs = m_perDofs2_u.size();
        m_perPdofs = m_perDofs2_p.size();

        m_dofs = m_blockAssemblerRef.numDofs() - m_blockAssemblerRef.getTarDim() * m_perUdofs - m_perPdofs;
        m_udofs = blockUdofs - m_perUdofs;
        m_pdofs = blockPdofs - m_perPdofs;
        m_pshift = m_blockAssemblerRef.getTarDim() * m_udofs;

        m_uMap = gsVector<index_t>::LinSpaced(blockUdofs, 0, blockUdofs - 1);
        m_uInvMap = m_uMap.transpose();
        for (std::vector<index_t>::iterator it = m_perDofs2_u.begin(); it != m_perDofs2_u.end(); ++it) {
            m_uInvMap.removeCol(m_uMap(*it));
            m_uMap.bottomRows(blockUdofs - (*it) - 1) -= gsVector<index_t>::Ones(blockUdofs - (*it) - 1);
        }
        m_uInvMap.transposeInPlace();
        for (std::vector<index_t>::iterator it = m_perDofs2_u.begin(); it != m_perDofs2_u.end(); ++it) {
            m_uMap(*it) = m_uMap(m_perDofs1_u[it - m_perDofs2_u.begin()]);
        }

        m_pMap = gsVector<index_t>::LinSpaced(blockPdofs, 0, blockPdofs - 1);
        m_pInvMap = m_pMap.transpose();
        for (std::vector<index_t>::iterator it = m_perDofs2_p.begin(); it != m_perDofs2_p.end(); ++it) {
            m_pInvMap.removeCol(m_pMap(*it));
            m_pMap.bottomRows(blockPdofs - (*it) - 1) -= gsVector<index_t>::Ones(blockPdofs - (*it) - 1);
        }
        m_pInvMap.transposeInPlace();
        for (std::vector<index_t>::iterator it = m_perDofs2_p.begin(); it != m_perDofs2_p.end(); ++it) {
            m_pMap(*it) = m_pMap(m_perDofs1_p[it - m_perDofs2_p.begin()]);
        }
    }

    // nonperiodic dof to periodic dof map
    inline index_t map(int i) const
    {
        if (i < m_blockAssemblerRef.getPshift())
            return (m_uMap(i % m_blockAssemblerRef.getUdofs()) + (i / m_blockAssemblerRef.getUdofs()) * m_udofs);
        else
        {
            return (m_pMap(i % m_blockAssemblerRef.getPshift()) + m_pshift);
        }
    }

    // periodic dof to nonperiodic dof map
    inline index_t invMap(int i) const
    {
        if (i < m_pshift)
            return (m_uInvMap(i % m_udofs) + (i / m_udofs) * m_blockAssemblerRef.getUdofs());
        else
            return (m_pInvMap(i % m_pshift) + m_blockAssemblerRef.getPshift());
    }

    // returns true if nonperiodic dof is eliminated 
    inline bool isEliminated(int i) const
    {
        if (i < m_blockAssemblerRef.getPshift())
            return (m_uIsEliminated(i % m_blockAssemblerRef.getUdofs()));
        else
            return (m_pIsEliminated(i % m_blockAssemblerRef.getPshift()));
    }

    // nonperiodic column coefficient vector to periodic coefficient vector
    void nonper2per_into(const gsMatrix<T> & nonperVector, gsMatrix<T> & perVector) const
    {
        perVector.setZero(m_dofs, 1);

        for (index_t row = 0; row < m_blockAssemblerRef.numDofs(); ++row)
            if (!isEliminated(row))
                perVector(map(row), 0) += nonperVector(row, 0);
            else {
                if (row < m_blockAssemblerRef.getPshift()) {
                    for (int t = 0; t < m_blockAssemblerRef.getTarDim(); t++)
                        perVector(m_uMap(row % m_blockAssemblerRef.getUdofs()) + t * m_udofs, 0) += m_transformMatrix(t, row / m_blockAssemblerRef.getUdofs()) * nonperVector(row, 0);
                }
                else
                    perVector(map(row), 0) += nonperVector(row, 0);
            }
    }

    // periodic column coefficient vector to nonperiodic coefficient vector
    void per2nonper_into(const gsMatrix<T> & perVector, gsMatrix<T> & nonperVector) const
    {
        nonperVector.setZero(m_blockAssemblerRef.numDofs(), 1);

        for (index_t row = 0; row < m_dofs; ++row)
            nonperVector(invMap(row), 0) += perVector(row, 0);

        for (int i = 0; i < m_perUdofs; i++) {
            gsMatrix<T> u1(m_blockAssemblerRef.getTarDim(), 1);
            for (int s = 0; s < m_blockAssemblerRef.getTarDim(); s++) {
                u1(s, 0) = perVector(m_uMap(m_perDofs2_u[i]) + s * m_udofs, 0);
            }
            gsMatrix<T> u2(m_blockAssemblerRef.getTarDim(), 1);
            u2 = m_invTransformMatrix * u1;
            for (int t = 0; t < m_blockAssemblerRef.getTarDim(); t++) {
                nonperVector(m_perDofs2_u[i] + t * m_blockAssemblerRef.getUdofs(), 0) = u2(t, 0);
            }
        }

        for (int i = 0; i < m_perPdofs; i++)
            nonperVector(m_perDofs2_p[i] + m_blockAssemblerRef.getPshift(), 0) = perVector(map(m_perDofs2_p[i] + m_blockAssemblerRef.getPshift()), 0);
    }

    void applyPeriodicConditionsOnVelocity_into(const gsSparseMatrix<T> & matrix, gsSparseMatrix<T> & resultMatrix, bool reserveMemory = true)
    {
        if (reserveMemory)
        {
            int pshiftNonper = m_blockAssemblerRef.getPshift();
            gsVector<int> nonZerosPerColumnVector;
            nonZerosPerColumnVector.setZero(m_pshift);

            for (int j = 0; j < pshiftNonper; j++)
                nonZerosPerColumnVector(map(j)) = matrix.col(j).nonZeros() + getPerUDofs();

            resultMatrix.resize(m_pshift, m_pshift);
            resultMatrix.reserve(nonZerosPerColumnVector);
        }

        short_t tarDim = m_blockAssemblerRef.getTarDim();
        int udofs = m_blockAssemblerRef.getUdofs();

        // MATRIX TRANSFORMATION      
        // map non-eliminated elements
        #pragma omp parallel for num_threads(m_blockAssemblerRef.getNumThreads())
        for (index_t col = 0; col < m_blockAssemblerRef.getPshift(); ++col)
            if (!(isEliminated(col)))
                for (typename gsSparseMatrix<T>::InnerIterator it(matrix, col); it; ++it)
                    if (!(isEliminated(it.row())))
                        resultMatrix.insert(map(it.row()), map(it.col())) = it.value();

        //col additions
        //U
        #pragma omp parallel for num_threads(m_blockAssemblerRef.getNumThreads())
        for (int i = 0; i < this->getPerUDofs(); i++)
            for (int s = 0; s < tarDim; s++)
                for (typename gsSparseMatrix<T>::InnerIterator it(matrix, m_perDofs2_u[i] + s * udofs); it; ++it)
                    if (!(isEliminated(it.row())))
                        for (int t = 0; t < tarDim; t++)
                            resultMatrix.coeffRef(map(it.row()), m_uMap(m_perDofs2_u[i]) + t * getUDofsPer()) += m_transformMatrix(t, s) * it.value();
                    //cross elements
                    else
                        for (int t = 0; t < tarDim; t++)
                            for (int u = 0; u < tarDim; u++)
                                resultMatrix.coeffRef(m_uMap(it.row() % udofs) + u * getUDofsPer(), m_uMap(m_perDofs2_u[i]) + t * getUDofsPer()) += m_transformMatrix(u, it.row() / udofs) * m_transformMatrix(t, s) * it.value();

        //row additions
        #pragma omp parallel for num_threads(m_blockAssemblerRef.getNumThreads())
        for (index_t col = 0; col < m_blockAssemblerRef.getPshift(); ++col)
            if (!(isEliminated(col)))
                for (typename gsSparseMatrix<T>::InnerIterator it(matrix, col); it; ++it)
                    if (isEliminated(it.row()))
                        for (int t = 0; t < tarDim; t++)
                            resultMatrix.coeffRef(m_uMap(it.row() % udofs) + t * getUDofsPer(), map(col)) += m_transformMatrix(t, it.row() / udofs) * it.value();

        resultMatrix.makeCompressed();
    }

    void applyPeriodicConditionsOnPressure_into(const gsSparseMatrix<T> & matrix, gsSparseMatrix<T> & resultMatrix, bool reserveMemory = true)
    {
        if (reserveMemory)
        {
            int pdofsNonper = m_blockAssemblerRef.getPdofs();
            gsVector<int> nonZerosPerColumnVector;
            nonZerosPerColumnVector.setZero(m_pdofs);

            for (int j = 0; j < pdofsNonper; j++)
                nonZerosPerColumnVector(m_pMap(j)) = matrix.col(j).nonZeros() + getPerPDofs();

            resultMatrix.resize(m_pdofs, m_pdofs);
            resultMatrix.reserve(nonZerosPerColumnVector);
        }

        for (index_t col = 0; col < matrix.cols(); ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(matrix, col); it; ++it)
                resultMatrix.coeffRef(m_pMap(it.row()), m_pMap(it.col())) += it.value();

        resultMatrix.makeCompressed();
    }

    gsField<T> constructSolution(const gsMatrix<T>& perSolVector, int unk, bool relative = false) const
    {
        gsMatrix<T> nonperSolVector;
        per2nonper_into(perSolVector, nonperSolVector);

        return m_blockAssemblerRef.constructSolution(nonperSolVector, unk, relative);
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& perSolVector, int unk, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, bool relative = false) const
    {
        gsMatrix<T> nonperSolVector;
        per2nonper_into(perSolVector, nonperSolVector);

        return m_blockAssemblerRef.constructSolution(nonperSolVector, unk, patchNumbers, sides, relative);
    }  

    gsField<T> constructSolutionCombined(const gsMatrix<T>& perSolVector, gsVector<size_t> relPatches) const
    {
        gsMatrix<T> nonperSolVector;
        per2nonper_into(perSolVector, nonperSolVector);

        return m_blockAssemblerRef.constructSolutionCombined(nonperSolVector, relPatches);
    }

    // computes flow rate through a side of a patch
    virtual T computeFlowRate(int patch, boxSide side, gsMatrix<T> perSolVector) const
    {
        gsMatrix<T> nonperSolVector;
        per2nonper_into(perSolVector, nonperSolVector);

        return m_blockAssemblerRef.computeFlowRate(patch, side, nonperSolVector);
    }

    virtual T computeDimensionlessWallDistance(gsMatrix<T> perSolVector, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        gsMatrix<T> nonperSolVector;
        per2nonper_into(perSolVector, nonperSolVector);

        return m_blockAssemblerRef.computeDimensionlessWallDistance(nonperSolVector, distancePatches, distanceSides, viscosity, reynoldsNumber, uFreeStream, maxYplus, npts, print, estimate);
    }

public:
    int numDofsPer() const { return m_dofs; }
    int getUDofsPer() const { return m_udofs; }
    int getPDofsPer() const { return m_pdofs; }
    int getPShiftPer() const { return m_pshift; }
    int getPerUDofs() const { return m_perUdofs; }
    int getPerPDofs() const { return m_perPdofs; }
    const gsMatrix<T>& getTransformMatrix() const { return m_transformMatrix; }
    const gsMatrix<T>& getInvTransformMatrix() const { return m_invTransformMatrix; }

protected:
    int m_dofs;
    int m_udofs;
    int m_pdofs;
    int m_pshift;

    gsMatrix<T> m_transformMatrix;
    gsMatrix<T> m_invTransformMatrix;

    int m_perUdofs, m_perPdofs;
    std::vector<index_t> m_perDofs1_u, m_perDofs2_u, m_perDofs1_p, m_perDofs2_p;
    gsMatrix<index_t> m_uMap, m_pMap;
    gsMatrix<index_t> m_uInvMap, m_pInvMap;
    gsMatrix<bool> m_uIsEliminated, m_pIsEliminated;

    bool m_bNonZerosExact;
    gsVector<int> m_nonZerosPerColumnVector;

    gsMatrix<T> m_perSolution;

    const uwbINSBlockAssembler<T>& m_blockAssemblerRef;

}; // class uwbINSAssemblerPeriodicBase

} // namespace gismo
