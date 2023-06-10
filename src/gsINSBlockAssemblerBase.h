/** @file gsINSBlockAssemblerBase.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSSolverParams.h>
#include <gsIncompressibleFlow/src/gsINSBlockVisitors.h>
//#include <gsIncompressibleFlow/src/gsINSdgBlockVisitors.h>

#include <gsCore/gsFuncData.h>

#ifdef _OPENMP 
#include <omp.h>
#endif

namespace gismo
{

/// @brief A base class for the assemblers of matrix blocks needed for incompressible flow problems.
/// @tparam T coefficient type
template<class T>
class gsINSBlockAssemblerBase
{

protected: // *** Class members ***

    int m_dofs;
    int m_numThreads;
    int m_tarDim;
    T m_viscosity;
    bool m_bUnsteady;
    gsINSSolverParams<T> m_params;
    std::vector<gsDofMapper> m_dofMappers;
    std::vector<gsMatrix<T> > m_ddof;
    gsMatrix<T> m_solution;
    std::vector< std::pair<index_t, index_t> > m_elementList, m_iFaceElementList;
    std::vector<boundaryInterface> m_iFaceList;

public:  // *** Constructor/destructor ***

    /// @brief Constructor of the object.
    /// @param[in] params a set of all needed parameters
    gsINSBlockAssemblerBase(const gsINSSolverParams<T>& params) : m_params(params)
    {
        initMembers();
    }

    virtual ~gsINSBlockAssemblerBase()
    {}

protected: // *** Member functions ***

    /// @brief Initialize the class members.
    virtual void initMembers();

    /**
     * @brief Create list of all elements for purpose of OMP parallelization. 
     * 
     * The first column of the resulting list contains patch indices and 
     * the second column contains element indices on the given patch.
    */
    void initElementList();


    /**
     * @brief Create list of all interface for purpose of OMP parallelization. 
     * 
     * The first column of the resulting list contains patch indices and 
     * the second column contains element indices on the given patch.
    */
    void initiFaceElementList();


    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  basisID     the index of the basis corresponding to \a unk (same as \a unk in this case)
    /// @param[in]  ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofs(const int unk, const int basisID, gsMatrix<T>& ddofVector);


    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using interpolation.
    /// @param[in]  unk          the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper       reference to the DOF mapper for \a unk
    /// @param[in]  mbasis       reference to the basis corresponding to \a unk
    /// @param[out] ddofVector   reference to the vector where computed coefficients will be stored
    void computeDirichletDofsIntpl(const int unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);


    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using L2-projection.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper      reference to the DOF mapper for \a unk
    /// @param[in]  mbasis      reference to the basis corresponding to \a unk
    /// @param[out] ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofsL2Proj(const int unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);


    /// @brief The main procedure for assembly of a block.
    /// @tparam     ElementBlockVisitor     type of the element visitor to be used for assembly (derived from gsINSBlockVisitor)
    /// @param[out] matrixBlock             the resulting matrix
    /// @param[out] rhs                     the resulting right-hand-side vector arising from eliminated Dirichlet DOFs
    /// @param[in]  nonZerosPerCol          number of nonzero elements per column for memory allocation
    template<class ElementBlockVisitor>
    void assembleBlock(gsSparseMatrix<T> & matrixBlock, gsMatrix<T> & rhs, const int nonZerosPerCol);
    
    /// @brief Initialize the element visitor according to its type.
    /// @param[in,out]  blockVisitor    a pointer to the block element visitor
    virtual void blockSpecificSettings(gsINSVisitorBase<T>* blockVisitor)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief The procedure for assembly of a right-hand side.
    /// @tparam         RhsVisitor  type of the element visitor to be used for assembly
    /// @param[in,out]  rhs         the resulting right-hand-side vector
    template<class RhsVisitor>
    void assembleRhs(gsMatrix<T> & rhs);

    /// @brief Initialize the RHS element visitor according to its type.
    /// @param[in,out]  blockVisitor    a pointer to the element visitor
    virtual void rhsSpecificSettings(gsINSVisitorBase<T>* blockVisitor)
    { GISMO_NO_IMPLEMENTATION }

    /* /// @brief The procedure for assembly of the block resulting from discontinuous Galerkin (DG) along the interfaces.
    /// @tparam     InterfaceBlockVisitor   type of the DG element visitor to be used for assembly (derived from gsINSBlockVisitorDG)
    /// @param[out] matrixBlock             the resulting matrix
    /// @param[out] rhs                     the resulting right-hand-side vector arising from eliminated Dirichlet DOFs
    /// @param[in]  nonZerosPerCol          number of nonzero elements per column for memory allocation */
    // template<class InterfaceBlockVisitor>
    // void assembleBlockDg(gsSparseMatrix<T> & matrixBlock, gsMatrix<T> & rhs, const int nonZerosPerCol);

public: // *** Member functions ***

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown (0 - velocity, 1 - pressure)
    void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0);

public: // *** Getters/setters ***

    /// @brief Returns the number of degrees of freedom (DOFs).
    int numDofs() const 
    { 
        GISMO_ASSERT(m_dofs > 0, "Something went wrong, number of DOFs is zero!");
        return m_dofs; 
    }

    /// @brief Returns the target dimension.
    int getTarDim() const { return m_tarDim; }

    /// @brief Returns the number of threads to be used by OpenMP.
    int getNumThreads() const { return m_numThreads; }

    /// @brief Set the number of threads to be used by OpenMP.
    void setNumThreads(const int numThreads)
    {
        #ifdef _OPENMP 
        const int maxThreads = omp_get_max_threads();
        if ((numThreads > 0) && (numThreads <= maxThreads))
            m_numThreads = numThreads;
        else
        {
            gsWarn << "The maximum number of threads ( " << maxThreads << " ) will be used.\n";
            m_numThreads = maxThreads;
        }
        #else
        m_numThreads = 1;
        #endif
    }

    /// @brief Returns true if the solved problem is unsteady.
    bool isUnsteady() const { return m_bUnsteady; }

    /// @brief Set the unsteady flag to the given value (true/false).
    void setUnsteady(bool unsteady = true) { m_bUnsteady = unsteady; }

    /// @brief Returns the viscosity value.
    T getViscosity() const { return m_viscosity; }

    /**
     * @brief Returns a const reference to the DOF mappers.
     *
     * The mapper for velocity is stored first, the mapper for pressure is second.
     */
    const std::vector<gsDofMapper>& getMappers() const { return m_dofMappers; }
    
    /**
     * @brief Returns a const reference to the vectors of coefficients at the Dirichlet boundaries.
     *
     * The vector of velocity coefficients is stored first, the vector of pressure coefficients is second.
     */
    const std::vector<gsMatrix<T> >& getDirichletDofs() const { return m_ddof; }
    
    /// @brief Returns a const reference to the current computed solution.
    const gsMatrix<T>& getSolution() const { return m_solution; }
    
    /// @brief Returns a const reference to the multipatch representing the computational domain.
    const gsMultiPatch<T>& getPatches() const { return m_params.getPde().patches(); }
    
    /**
     * @brief Returns a reference to the discretization bases.
     *
     * The velocity basis is stored first, the  pressure basis is second.
     * 
     * There is also a const version returning a const reference.
     */
    std::vector< gsMultiBasis<T> >& getBases() { return m_params.getBases(); }
    const std::vector< gsMultiBasis<T> >& getBases() const { return m_params.getBases(); }
    
    /// @brief Returns a const reference to the boundary conditions.
    const gsBoundaryConditions<T>& getBCs() const { return m_params.getBCs(); }
    
    /// @brief Returns a pointer to the right-hand-side function.
    const gsFunction<T>* getRhsFcn() const { return m_params.getPde().rhs(); }

    /// @brief Returns the assembler options.
    gsAssemblerOptions  getAssemblerOptions() const { return m_params.assemblerOptions(); }
    
    /// @brief Returns the INS solver option list.
    gsOptionList    options() const { return m_params.options(); }

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSBlockAssemblerBase.hpp)
#endif
