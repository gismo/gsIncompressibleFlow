/** @file gsFlowAssemblerBase.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowVisitors.h>

namespace gismo
{

/// @brief  A           base class for all assemblers in gsIncompressibleFlow.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsFlowAssemblerBase
{

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    index_t m_dofs;
    short_t m_tarDim;
    bool m_isInitialized;
    bool m_isBaseReady;
    bool m_isSystemReady;
    std::vector<gsMatrix<T> > m_ddof;
    gsMatrix<T> m_solution;


public: // *** Constructor/destructor ***

    gsFlowAssemblerBase(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    m_paramsPtr(paramsPtr)
    { }

    virtual ~gsFlowAssemblerBase()
    { }


protected: // *** Member functions ***

    /// @brief Initialize the class members.
    void initMembers();

    /// @brief Update sizes of members (when DOF numbers change after constructing the assembler).
    virtual void updateSizes()
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  basisID     the index of the basis corresponding to \a unk (same as \a unk in this case)
    /// @param[in]  ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofs(const index_t unk, const index_t basisID, gsMatrix<T>& ddofVector);

    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using interpolation.
    /// @param[in]  unk          the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper       reference to the DOF mapper for \a unk
    /// @param[in]  mbasis       reference to the basis corresponding to \a unk
    /// @param[out] ddofVector   reference to the vector where computed coefficients will be stored
    void computeDirichletDofsIntpl(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);

    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using L2-projection.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper      reference to the DOF mapper for \a unk
    /// @param[in]  mbasis      reference to the basis corresponding to \a unk
    /// @param[out] ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofsL2Proj(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);

    /// @brief Assemble a matrix block.
    /// @param[in]  visitor     visitor for the required block
    /// @param[in]  testBasisID ID of the test basis
    /// @param[out] block       the resulting matrix block
    /// @param[out] blockRhs    right-hand side for the matrix block (arising from eliminated Dirichlet DOFs)
    /// @param[in]  compressMat call makeCompressed() at the end
    template<class ElementVisitor>
    void assembleBlock(ElementVisitor& visitor, index_t testBasisID, gsSparseMatrix<T, MatOrder>& block, gsMatrix<T>& blockRhs, bool compressMat = true);

    /// @brief Assemble the right-hand side.
    /// @param[in]  visitor     visitor for the right-hand side
    /// @param[in]  testBasisID ID of the test basis
    /// @param[out] rhs         the resulting right-hand side vector
    template<class ElementVisitorRhs>
    void assembleRhs(ElementVisitorRhs& visitor, index_t testBasisID, gsMatrix<T>& rhs);

    /// @brief Assemble the linear part of the problem.
    virtual void assembleLinearPart()
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Assemble the nonlinear part of the problem.
    virtual void assembleNonlinearPart()
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Assemble all that needs to be updated in each nonlinear iteration.
    virtual void updateAssembly();


    /// @brief Update current solution field stored in the assembler.
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
    {GISMO_NO_IMPLEMENTATION}


public: // *** Member functions ***

    /// @brief Initialize the assembler.
    virtual void initialize();

    /// @brief Update the assembler in new nonlinear iteration.
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution
    virtual void update(const gsMatrix<T> & solVector, bool updateSol = true);

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector    the solution vector obtained from the linear system
    /// @param[out] result       the resulting solution as a gsMultiPatch object
    /// @param[in]  unk          the considered unknown
    /// @param[in]  customSwitch a switch to be used for any purpose by derived classes
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, index_t unk, bool customSwitch = false) const
    {GISMO_NO_IMPLEMENTATION}


public: // *** Getters/setters ***

    /// @brief Returns the number of degrees of freedom (DOFs).
    index_t numDofs() const 
    { 
        GISMO_ASSERT(m_dofs > 0, "Something went wrong, number of DOFs is zero!");
        return m_dofs; 
    }

    /// @brief Returns the target dimension.
    short_t getTarDim() const { return m_tarDim; }

    /// @brief Returns true if the assembler has been initialized.
    bool isInitialized() { return m_isInitialized; }

    /**
     * @brief Returns a const reference to the vectors of coefficients at the Dirichlet boundaries.
     *
     * In the case of velocity and pressure, the vector of velocity coefficients is stored first, the vector of pressure coefficients is second.
     */
    const std::vector<gsMatrix<T> >& getDirichletDofs() const { return m_ddof; }

    /// @brief Returns a const reference to the current computed solution.
    const gsMatrix<T>& getSolution() const { return m_solution; }

    /// @brief Returns a const reference to the multipatch representing the computational domain.
    const gsMultiPatch<T>& getPatches() const { return m_paramsPtr->getPde().patches(); }

    /**
     * @brief Returns a reference to the discretization bases.
     *
     * Order of bases: velocity, pressure, (turb. model quantities)
     * 
     * There is also a const version returning a const reference.
     */
    virtual std::vector< gsMultiBasis<T> >& getBases() { return m_paramsPtr->getBases(); }
    virtual const std::vector< gsMultiBasis<T> >& getBases() const { return m_paramsPtr->getBases(); }

    /**
     * @brief Returns a reference to the discretization bases for variable \a unk.
     *
     * Order of bases: velocity, pressure, (turb. model quantities)
     * 
     * @param[in] unk unknown index
     * 
     * There is also a const version returning a const reference.
     */
    virtual gsMultiBasis<T>& getBasis(index_t unk) { return m_paramsPtr->getBasis(unk); }
    virtual const gsMultiBasis<T>& getBasis(index_t unk) const { return m_paramsPtr->getBasis(unk); }

    /**
     * @brief Returns a reference to the DOF mappers.
     *
     * Order of mappers: velocity, pressure, (turb. model quantities)
     * 
     * There is also a const version returning a const reference.
     */
    virtual std::vector< gsDofMapper >& getMappers() { return m_paramsPtr->getMappers(); }
    virtual const std::vector< gsDofMapper >& getMappers() const { return m_paramsPtr->getMappers(); }

    /**
     * @brief Returns a reference to the DOF mapper for variable \a unk.
     *
     * Order of mappers: velocity, pressure, (turb. model quantities)
     * 
     * @param[in] unk unknown index
     * 
     * There is also a const version returning a const reference.
     */
    virtual gsDofMapper& getMapper(index_t unk) { return m_paramsPtr->getMapper(unk); }
    virtual const gsDofMapper& getMapper(index_t unk) const { return m_paramsPtr->getMapper(unk); }

    /// @brief Returns a const reference to the boundary conditions.
    virtual const gsBoundaryConditions<T>& getBCs() const { return m_paramsPtr->getBCs(); }

    /// @brief Returns a pointer to the right-hand-side function.
    const gsFunction<T>* getRhsFcn() const { return m_paramsPtr->getPde().rhs(); }

    /// @brief Returns the assembler options.
    gsAssemblerOptions getAssemblerOptions() const { return m_paramsPtr->assemblerOptions(); }
    
    /// @brief Returns the flow solver option list.
    gsOptionList options() const { return m_paramsPtr->options(); }

    /// @brief Returns the assembled matrix.
    virtual const gsSparseMatrix<T, MatOrder>& matrix() const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the assembled matrix for unknown with index \a unk (e.g. from two-equation turbulence models).
    /// @param[in] unk index of the unknown
    virtual const gsSparseMatrix<T, MatOrder>& matrix(index_t unk) const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the mass matrix for unknown with index \a unk. There is also a const version.
    /// @param[in] unkID index of the unknown
    virtual gsSparseMatrix<T, MatOrder>& getMassMatrix(index_t unkID)
    {GISMO_NO_IMPLEMENTATION}
    
    virtual const gsSparseMatrix<T, MatOrder>& getMassMatrix(index_t unkID) const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the assembled right-hand side.
    virtual const gsMatrix<T>& rhs() const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the assembled right-hand side for unknown with index \a unk (e.g. from two-equation turbulence models).
    /// @param[in] unk index of the unknown
    virtual const gsMatrix<T>& rhs(index_t unk) const
    {GISMO_NO_IMPLEMENTATION}

};

// ---------------------------------------------------------------------------
// definitions of template member functions

template<class T, int MatOrder>
template<class ElementVisitor>
void gsFlowAssemblerBase<T, MatOrder>::assembleBlock(ElementVisitor& visitor, index_t testBasisID, gsSparseMatrix<T, MatOrder>& block, gsMatrix<T>& blockRhs, bool compressMat)
{
    if (m_paramsPtr->options().getString("assemb.loop") == "RbR")
    {
        for(size_t p = 0; p < getPatches().nPatches(); p++)
        {
            visitor.initOnPatch(p);
            index_t nBases = m_paramsPtr->getBasis(testBasisID).piece(p).size();

            for(index_t i = 0; i < nBases; i++)
            {
                visitor.evaluate(i);
                visitor.assemble();
                visitor.localToGlobal(m_ddof, block, blockRhs);
            }
        }
    }
    else        
    {
        if (m_paramsPtr->options().getString("assemb.loop") != "EbE")
            gsWarn << "Unknown matrix formation method, using EbE (element by element)!\n";

        #pragma omp parallel
        { 
            ElementVisitor
            #ifdef _OPENMP
            // Create thread-private visitor
            visitor_(visitor);
            const int threadId = omp_get_thread_num();
            const int numThreads  = omp_get_num_threads();
            #else
            &visitor_ = visitor;
            #endif
        
            // iteration over all elements in all patches
            typename gsBasis<T>::domainIter domIt = m_paramsPtr->getBasis(testBasisID).domain()->beginAll();
            typename gsBasis<T>::domainIter domItEnd = m_paramsPtr->getBasis(testBasisID).domain()->endAll();

            index_t patchID = -1;

            #ifdef _OPENMP
            domIt += threadId;
            for (; domIt < domItEnd; domIt+=(numThreads) )
            #else
            for (; domIt < domItEnd; ++domIt )
            #endif
            {
                index_t p = domIt.patch();
                if (p != patchID)
                {
                    patchID = p;
                    visitor_.initOnPatch(patchID);
                }

                visitor_.evaluate(domIt.get());
                visitor_.assemble();

                #pragma omp critical(localToGlobal) // only one thread at a time
                visitor_.localToGlobal(m_ddof, block, blockRhs);
            }
        }
    }   
}


template<class T, int MatOrder>
template<class ElementVisitorRhs>
void gsFlowAssemblerBase<T, MatOrder>::assembleRhs(ElementVisitorRhs& visitor, index_t testBasisID, gsMatrix<T>& rhs)
{
    if (m_paramsPtr->options().getString("assemb.loop") == "RbR")
    {
        for(size_t p = 0; p < getPatches().nPatches(); p++)
        {
            visitor.initOnPatch(p);
            index_t nBases = m_paramsPtr->getBasis(testBasisID).piece(p).size();

            for(index_t i = 0; i < nBases; i++)
            {
                visitor.evaluate(i);
                visitor.assemble();
                visitor.localToGlobal(rhs);
            }
        }
    }
    else        
    {
        if (m_paramsPtr->options().getString("assemb.loop") != "EbE")
            gsWarn << "Unknown matrix formation method, using EbE (element by element)!\n";

        #pragma omp parallel
        { 
            ElementVisitorRhs
            #ifdef _OPENMP
            // Create thread-private visitor
            visitor_(visitor);
            const int threadId = omp_get_thread_num();
            const int numThreads  = omp_get_num_threads();
            #else
            &visitor_ = visitor;
            #endif
        
            // iteration over all elements in all patches
            typename gsBasis<T>::domainIter domIt = m_paramsPtr->getBasis(testBasisID).domain()->beginAll();
            typename gsBasis<T>::domainIter domItEnd = m_paramsPtr->getBasis(testBasisID).domain()->endAll();

            index_t patchID = -1;

            #ifdef _OPENMP
            domIt += threadId;
            for (; domIt < domItEnd; domIt+=(numThreads) )
            #else
            for (; domIt < domItEnd; ++domIt )
            #endif
            {
                index_t p = domIt.patch();
                if (p != patchID)
                {
                    patchID = p;
                    visitor_.initOnPatch(patchID);
                }

                visitor_.evaluate(domIt.get());
                visitor_.assemble();

                #pragma omp critical(localToGlobal) // only one thread at a time
                visitor_.localToGlobal(rhs);
            }
        }
    }   
}


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowAssemblerBase.hpp)
#endif