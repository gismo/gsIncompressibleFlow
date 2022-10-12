/** @file uwbADRAssemblerBase.h

Author(s): E. Turnerova
*/

#pragma once

#ifdef _OPENMP 
#include <omp.h>
#endif

#include "uwbADRBlockAssembler.h"

namespace gismo
{

template<class T>
class uwbADRSolverParams;

template<class T>
class uwbADRAssemblerSteady
{
public:
    uwbADRAssemblerSteady(uwbADRSolverParams<T>& params, int numVar): m_blockAssembler(params, numVar)
    {
        initMembers();
    }

    virtual ~uwbADRAssemblerSteady()
    {
    }

protected:
    void initMembers()
    {
        int numVar = getNumVar();
        int varDofs = numVarDofs();
        int dofs = numDofs();

        m_matrix.resize(varDofs, dofs);
        m_rhs.setZero(varDofs, numVar);
        m_baseMatrix.resize(varDofs, dofs);
        m_baseRhs.setZero(varDofs, numVar);

        m_solution.setZero(varDofs, numVar);

        m_bInitialized = false;
        m_bSystemReady = false;
    }

public:
    virtual void initialize()
    {
        initSteadyAssembly();

        fillBase();

        m_bInitialized = true;
    }

    void setInitialCondition(const gsMatrix<T> & solVector)
    {
        m_solution = solVector;
        m_blockAssembler.setSolution(solVector);
    }

    virtual void updateLinearSteady()
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_rhs = m_baseRhs;
    }

    virtual void updatePicardSteady(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        updatePicardAssembly(solVector);

        fillSystem();
        m_rhs = m_baseRhs;
    }

protected:
    virtual void initSteadyAssembly() { GISMO_NO_IMPLEMENTATION }

    virtual void updatePicardAssembly(const gsMatrix<T> & solVector) { GISMO_NO_IMPLEMENTATION }

    void checkReady() const { GISMO_ASSERT(m_bSystemReady, "System not ready, update() must be called first"); }

public:

    virtual void updateLinear(const gsMatrix<T> & solVector) { GISMO_NO_IMPLEMENTATION }

    virtual void fillBase() { GISMO_NO_IMPLEMENTATION }

    virtual void fillSystem() { GISMO_NO_IMPLEMENTATION }

    const gsSparseMatrix<T> matrix(index_t var = 0) const
    {
        GISMO_ASSERT(var < getNumVar(), "Matrix index higher than number of variables.");
        checkReady();

        return m_matrix.middleCols(var * this->numVarDofs(), this->numVarDofs());
    }

    const gsMatrix<T> rhs(index_t var = 0) const
    {
        GISMO_ASSERT(var < getNumVar(), "Rhs index higher than number of variables.");
        checkReady();

        return m_rhs.col(var);
    }

public:
    bool isInitialized() { return m_bInitialized; }
    const uwbADRBlockAssembler<T>& getBlockAssembler() const { return m_blockAssembler; }
    uwbADRBlockAssembler<T>& getBlockAssembler() { return m_blockAssembler; }
    int getNumVar() const { return m_blockAssembler.getNumVar(); }
    int numVarDofs() const { return m_blockAssembler.numVarDofs(); }
    int numDofs() const { return m_blockAssembler.numDofs(); }

    const gsMatrix<T>&  getSolution() const { return m_solution; }
    gsField<T> constructSolution(const gsMatrix<T>& solVector) const { return m_blockAssembler.constructSolution(solVector); }

protected:

    uwbADRBlockAssembler<T> m_blockAssembler;

    bool m_bInitialized;
    bool m_bSystemReady;

    // base system matrix, which does not change at all
    gsSparseMatrix<T> m_baseMatrix;
    gsMatrix<T> m_baseRhs;

    // complete system matrix, nBaseMatrix + the part changing with Picard iteration
    gsSparseMatrix<T> m_matrix;
    gsMatrix<T> m_rhs;

    gsMatrix<T> m_solution; // turbulent quantities solution

}; // class uwbADRAssemblerBase



//====================================================================================================================

template<class T>
class uwbADRAssemblerBaseUnsteady : public uwbADRAssemblerSteady<T>
{
public:
    typedef uwbADRAssemblerSteady<T> Base;

public:
    uwbADRAssemblerBaseUnsteady(uwbADRSolverParams<T>& params, int numVar) : Base(params, numVar)
    {
        m_timeStepSize = params.settings().get(constantsADR::timeStep);
    }

    virtual ~uwbADRAssemblerBaseUnsteady() { }

public:
    virtual void initialize()
    {
        initAssembly();

        this->fillBase();

        m_bInitialized = true;
    }

    virtual void updateLinear(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_solution = solVector;

        fillRhs();
    }

    virtual void update(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_solution = solVector;

        updateAssembly();

        this->fillSystem();
        fillRhs();
    }

    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        this->updatePicardAssembly(solVector);

        this->fillSystem();
        fillRhs();
    }

protected:
    virtual void initAssembly() { GISMO_NO_IMPLEMENTATION }
    virtual void fillRhs() { GISMO_NO_IMPLEMENTATION }
    virtual void updateAssembly() { GISMO_NO_IMPLEMENTATION }

public:
    const T getTimeStepSize() const { return m_timeStepSize; }
    
    void changeTimeStepSize(const T timeStepSize)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_timeStepSize = timeStepSize;

        m_baseMatrix.resize(numVarDofs(), numDofs());
        this->fillBase();

        m_bSystemReady = false;
    }

protected:

    T m_timeStepSize;

    // members from uwbADRAssemblerSteady
    using Base::m_bInitialized;
    using Base::m_baseMatrix;
    using Base::m_bSystemReady;
    using Base::m_blockAssembler;
    using Base::m_solution;
    using Base::m_baseRhs;

public:
    // functions from uwbADRAssemblerSteady
    using Base::numDofs;
    using Base::numVarDofs;
    using Base::getNumVar;

}; // class uwbADRAssemblerBaseUnsteady

} //namespace gismo
