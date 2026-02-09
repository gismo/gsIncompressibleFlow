/** @file gsINSVisitors.h

    @brief Incompressible Navier-Stokes visitors.
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsTMModels.h> 

namespace gismo
{

// ===================================================================================================================
// VELOCITY-VELOCITY VISITORS

/**
 * @brief Base visitor for the velocity-velocity part of the Navier-Stokes system.
 * 
 * Assumes identical matrix for all velocity components,
 * i.e. only one local matrix is evaluated and 
 * either global matrix for only one velocity component is assembled 
 * or the local matrix is distributed into all diagonal blocks of the global velocity matrix.
 * 
 * @tparam T           real number type
 * @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
 * @ingroup IncompressibleFlow
 */
template <class T, int MatOrder>
class gsINSVisitorUU : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_localMat;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;

public: // *** Constructor/destructor ***

    gsINSVisitorUU() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUU(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }
        

protected: // *** Member functions ***

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_trialUnkID = 0;   // velocity
    }

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);
};

// ===================================================================================================================

/// @brief Visitor for the linear terms in the velocity-velocity part of the Navier-Stokes system.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorUUlin : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUlin() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUUlin(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_Diffusion<T>(m_paramsPtr->getPde().viscosity()) );

        // ... other terms, e.g. from stabilizations
    }

};

// ===================================================================================================================

/// @brief Visitor for the non-linear terms in the velocity-velocity part of the Navier-Stokes system.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorUUnonlin : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


public:
     std::vector< gsMatrix<T> > m_TurbulentViscosityGrads;

protected: // *** Base class members ***

    using Base::m_solField;
    using Base::m_viscosity;
    using Base::m_TMModelPtr;
    using Base::m_TurbulentViscosityVals;
    using Base::m_tauS;
    using Base::m_USolVals;
    using Base::m_velocities;

    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_terms;
    using Base::m_quNodes;   

public: // *** Constructor/destructor ***

    gsINSVisitorUUnonlin() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUUnonlin(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { 
        m_viscosity = m_paramsPtr->getPde().viscosity();
    }

    /// @brief Copy constructor.
    gsINSVisitorUUnonlin(gsINSVisitorUUnonlin<T, MatOrder> const & other) :
    gsINSVisitorUUnonlin()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();

        //setting current velocity field in the terms
        this->setCurrentSolution(m_solField);
    }

protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new typename gsFlowVisitor<T, MatOrder>::ConvectionTerm() );

        // ... other terms, e.g. from stabilizations
        if ((m_paramsPtr->options().getSwitch("TCSD_NS")) || (m_paramsPtr->options().getSwitch("TCSD_RANS")) || (m_paramsPtr->options().getSwitch("SUPG_NS")) || (m_paramsPtr->options().getSwitch("SUPG_RANS")))
        {
            m_terms.push_back( new gsFlowTerm_TCSDStabilization_advection<T>() );
        }

        if ((m_paramsPtr->options().getSwitch("SUPG_NS")) || (m_paramsPtr->options().getSwitch("SUPG_RANS")))
        {
            m_terms.push_back( new gsFlowTerm_SUPGStabilization_diffusion<T>(m_viscosity) );
        }
    }

public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

public: // *** Getters/setters ***

    void setTurbulenceModel(typename gsTMModelData<T>::Ptr TMModelPtr) { m_TMModelPtr = TMModelPtr;}

};

// ===================================================================================================================

/// @brief Visitor for the velocity mass matrix.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorUUmass : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUmass() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUUmass(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_ValVal<T>() );
    }

};

// ===================================================================================================================

/// @brief Visitor for the terms arising from (omega x r) in rotating frame of reference.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorUUrotation : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;

protected: // *** Class members ***

    real_t m_omega;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_patchID;
    using Base::m_localMat;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;


public: // *** Constructor/destructor ***

gsINSVisitorUUrotation() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUUrotation(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { 
        m_omega = m_paramsPtr->options().getReal("omega");
    }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        if (m_omega != 0)
            m_terms.push_back( new gsFlowTerm_ValVal<T>() );
    }

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

/// @brief Visitor for the time discretization term.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorUUtimeDiscr : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;

public: // *** Constructor/destructor ***

    gsINSVisitorUUtimeDiscr() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUUtimeDiscr(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );
    }

};

// ===================================================================================================================

/// @brief Visitor for the time discretization term.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorUU_TCSD_time : public gsINSVisitorUUnonlin<T, MatOrder>
{

public:
    typedef gsINSVisitorUUnonlin<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_solField;
    using Base::m_viscosity;
    using Base::m_TMModelPtr;
    using Base::m_TurbulentViscosityVals;
    using Base::m_tauS;
    using Base::m_USolVals;
    using Base::m_velocities;

    using Base::m_paramsPtr;
    using Base::m_terms;
    using Base::m_quNodes;

public: // *** Constructor/destructor ***

    gsINSVisitorUU_TCSD_time() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUU_TCSD_time(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }

    /// @brief Copy constructor.
    gsINSVisitorUU_TCSD_time(gsINSVisitorUU_TCSD_time<T, MatOrder> const & other) :
    gsINSVisitorUU_TCSD_time()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();

        //setting current velocity field in the terms
        this->setCurrentSolution(m_solField);

    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TCSDStabilization_time<T>(m_paramsPtr->options().getReal("timeStep")) );
    }

public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

public: // *** Getters/setters ***

    void setTurbulenceModel(typename gsTMModelData<T>::Ptr TMModelPtr) { m_TMModelPtr = TMModelPtr;}

};

// ===================================================================================================================
// ===================================================================================================================

// VELOCITY-PRESSURE VISITORS

/**
 * @brief Visitor for the velocity-pressure part of the Navier-Stokes system.
 * 
 * This visitor assembles for the terms with velocity test functions and pressure trial functions.
 * 
 * @tparam T            real number type
 * @tparam MatOrder     sparse matrix storage order (ColMajor/RowMajor)
 * 
 * @ingroup IncompressibleFlow
 */
template <class T, int MatOrder>
class gsINSVisitorPU : public gsFlowVisitorVectorValued<T, MatOrder>  // PU: trial, test
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;

public: // *** Constructor/destructor ***

    gsINSVisitorPU() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorPU(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTerm_PvalUdiv<T>() );
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_trialUnkID = 1;   // pressure
    }

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);
};

// ===================================================================================================================

/**
 * @brief Visitor for the velocity-pressure part of the Navier-Stokes system.
 * 
 * This visitor assembles for the terms with velocity test functions and pressure trial functions.
 * It also adds the right-hand side arising from the pressure-velocity part (from velocity Dirichlet boundary conditions) into the global right-hand side.
 * This is useful in the case, when the pressure-velocity part of the linear system is not assembled, since it is equal to the (negative) transpose of the velocity-pressure part.
 * 
 * @tparam T            real number type
 * @tparam MatOrder     sparse matrix storage order (ColMajor/RowMajor)
 * 
 * @ingroup IncompressibleFlow
 */
template <class T, int MatOrder>
class gsINSVisitorPU_withUPrhs : public gsINSVisitorPU<T, MatOrder>  // PU: trial, test
{

public:
    typedef gsINSVisitorPU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;

public: // *** Constructor/destructor ***

    gsINSVisitorPU_withUPrhs() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorPU_withUPrhs(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================
/**
 * @brief Visitor for the velocity-pressure part of the RANS equations with SUPG stabilization.
 * 
 * This visitor assembles for the terms with velocity test functions and pressure trial functions.
 * 
 * @tparam T            real number type
 * @tparam MatOrder     sparse matrix storage order (ColMajor/RowMajor)
 * 
 * @ingroup IncompressibleFlow
 */
/*template <class T, int MatOrder>
class gsINSVisitorPU_SUPG_presssure : public gsINSVisitorPU<T, MatOrder>
{

public:
    typedef gsINSVisitorPU<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_solField;
    using Base::m_viscosity;
    using Base::m_TMModelPtr;
    using Base::m_TurbulentViscosityVals;
    using Base::m_tauS;
    using Base::m_USolVals;
    using Base::m_velocities;

    using Base::m_paramsPtr;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_terms;
    using Base::m_quNodes;
    
public: // *** Constructor/destructor ***

    gsINSVisitorPU_SUPG_presssure() {}

    gsINSVisitorPU_SUPG_presssure(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    {   }

    /// @brief Copy constructor.
    gsINSVisitorPU_SUPG_presssure(gsINSVisitorPU_SUPG_presssure<T, MatOrder> const & other) :
    gsINSVisitorPU_SUPG_presssure()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();

        //setting current velocity field in the terms
        this->setCurrentSolution(m_solField);

    }

protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_SUPGStabilization_pressure<T>() );
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_trialUnkID = 1;   // pressure
    }

public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

public: // *** Getters/setters ***

    void setTurbulenceModel(typename gsTMModelData<T>::Ptr TMModelPtr) { m_TMModelPtr = TMModelPtr;} 

};
*/

// ===================================================================================================================
// ===================================================================================================================

// PRESSURE-VELOCITY VISITORS

/**
 * @brief Base visitor for the pressure-velocity part of the Navier-Stokes system.
 * 
 * These are visitors for the terms with pressure test functions and velocity trial functions.
 * 
 * @tparam T            real number type
 * @tparam MatOrder     sparse matrix storage order (ColMajor/RowMajor)
 * 
 * @ingroup IncompressibleFlow
 */
template <class T, int MatOrder>
class gsINSVisitorUP : public gsFlowVisitorVectorValued<T, MatOrder>  // UP: trial, test
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;

public: // *** Constructor/destructor ***

    gsINSVisitorUP() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorUP(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTerm_UdivPval<T>() );
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_trialUnkID = 0;   // velocity
    }

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);
};

// ===================================================================================================================
// ===================================================================================================================

// PRESSURE-PRESSURE VISITORS

/// @brief Base visitor for the pressure-pressure part of the Navier-Stokes system.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorPP : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_localMat;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;

public: // *** Constructor/destructor ***

    gsINSVisitorPP() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorPP(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_trialUnkID = 1;   // pressure
    }

protected: // *** Member functions ***

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
    {
        localToGlobal_nonper(eliminatedDofs, globalMat, globalRhs);
    }
};

// ===================================================================================================================

// template <class T, int MatOrder>
// class gsINSVisitorPPlin : public gsINSVisitorPP<T, MatOrder>
// {

// public:
//     typedef gsINSVisitorPP<T, MatOrder> Base;


// protected: // *** Base class members ***

//     using Base::m_paramsPtr;
//     using Base::m_terms;


// public: // *** Constructor/destructor ***

//     gsINSVisitorPPlin() {}

//     gsINSVisitorPPlin(typename gsFlowSolverParams<T>::Ptr paramsPtr):
//     Base(paramsPtr)
//     { }


// protected: // *** Member functions ***

//     virtual void defineTerms()
//     {
//         // no default pressure-pressure terms in the Navier-Stokes eqns
//         // optionally stabilization for inf-sup unstable discretizations
//     }

// };

// // ===================================================================================================================

// template <class T, int MatOrder>
// class gsINSVisitorPPnonlin : public gsINSVisitorPP<T, MatOrder>
// {

// public:
//     typedef gsINSVisitorPP<T, MatOrder> Base;


// protected: // *** Base class members ***

//     using Base::m_paramsPtr;
//     using Base::m_terms;


// public: // *** Constructor/destructor ***

//     gsINSVisitorPPnonlin() {}

//     gsINSVisitorPPnonlin(typename gsFlowSolverParams<T>::Ptr paramsPtr):
//     Base(paramsPtr)
//     { }


// protected: // *** Member functions ***

//     virtual void defineTerms()
//     {
//         // no default pressure-pressure terms in the Navier-Stokes eqns
//         // optionally stabilization for inf-sup unstable discretizations
//     }

// };

// ===================================================================================================================

/// @brief Visitor for the pressure mass matrix.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorPPmass : public gsINSVisitorPP<T, MatOrder>
{

public:
    typedef gsINSVisitorPP<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPmass() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorPPmass(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_ValVal<T>() );
    }

};

// ===================================================================================================================

/**
 * @brief Visitor for the pressure Laplacian operator.
 * 
 * This operator is needed for definition of the PCD preconditioner.
 * 
 * @tparam T            real number type
 * @tparam MatOrder     sparse matrix storage order (ColMajor/RowMajor)
 * 
 * @ingroup IncompressibleFlow
 */
template <class T, int MatOrder>
class gsINSVisitorPPlaplace : public gsINSVisitorPP<T, MatOrder>
{

public:
    typedef gsINSVisitorPP<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPlaplace() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorPPlaplace(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_GradGrad<T>() );
    }

};

// ===================================================================================================================

/**
 * @brief Visitor for the pressure convection operator.
 * 
 * This operator is needed for definition of the PCD preconditioner.
 * 
 * @tparam T            real number type
 * @tparam MatOrder     sparse matrix storage order (ColMajor/RowMajor)
 * 
 * @ingroup IncompressibleFlow
 */
template <class T, int MatOrder>
class gsINSVisitorPPconvection : public gsINSVisitorPP<T, MatOrder>
{

public:
    typedef gsINSVisitorPP<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPconvection() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorPPconvection(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTerm_UsolGradVal<T>() );
    }

};

// ===================================================================================================================
/**
 * @brief Visitor for the pressure-pressure stabilization term.
 * 
 * This visitor assembles for the residual-based stabilization term in the continuity equation.
 * 
 * @tparam T            real number type
 * @tparam MatOrder     sparse matrix storage order (ColMajor/RowMajor)
 * 
 * @ingroup IncompressibleFlow
 */
/*template <class T, int MatOrder>
class gsINSVisitorPP_ResidualStabilization_continuity : public gsINSVisitorPP<T, MatOrder>
{

public:
    typedef gsINSVisitorPP<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_solField;
    using Base::m_viscosity;
    using Base::m_TMModelPtr;
    using Base::m_TurbulentViscosityVals;
    using Base::m_tauS;
    using Base::m_USolVals;
    using Base::m_velocities;

    using Base::m_paramsPtr;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_terms;
    using Base::m_quNodes;
    
public: // *** Constructor/destructor ***

    gsINSVisitorPP_ResidualStabilization_continuity() {}

    gsINSVisitorPP_ResidualStabilization_continuity(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    {   }

    /// @brief Copy constructor.
    gsINSVisitorPP_ResidualStabilization_continuity(gsINSVisitorPP_ResidualStabilization_continuity<T, MatOrder> const & other) :
    gsINSVisitorPP_ResidualStabilization_continuity()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();

        //setting current velocity field in the terms
        this->setCurrentSolution(m_solField);

    }

protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_ResidualStabilization_continuity<T>() );
    }

public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

public: // *** Getters/setters ***

    void setTurbulenceModel(typename gsTMModelData<T>::Ptr TMModelPtr) { m_TMModelPtr = TMModelPtr;}

};
*/

// ===================================================================================================================

// template <class T, int MatOrder>
// class gsINSVisitorPP_PCDrobinBC : public gsINSVisitorPP<T, MatOrder>
// {

// public:
//     typedef gsINSVisitorPP<T, MatOrder> Base;

// protected: // *** Base class members ***

//     using Base::m_paramsPtr;
//     using Base::m_terms;


// public: // *** Constructor/destructor ***

//     gsINSVisitorPP_PCDrobinBC() {}

//     gsINSVisitorPP_PCDrobinBC(typename gsFlowSolverParams<T>::Ptr paramsPtr):
//     Base(paramsPtr)
//     { }


// protected: // *** Member functions ***

//     virtual void defineTerms()
//     {
//         // TODO
//     }

// };

// ===================================================================================================================
// ===================================================================================================================

// RHS VISITORS

/// @brief Visitor for right-hand side of the momentum equations (force function).
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorRhsU : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_localMat;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_trialFunData;
    using Base::m_periodicTransformMat;


public: // *** Constructor/destructor ***

    gsINSVisitorRhsU() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorRhsU(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr), m_pRhsFun(paramsPtr->getPde().force())
    {
        GISMO_ASSERT(m_pRhsFun->targetDim() == m_paramsPtr->getPde().dim(), "Wrong RHS function passed into gsINSVisitorRhsU.");
    }
        

protected: // *** Member functions ***

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_trialUnkID = 0;    // velocity (not needed here)
    }

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_rhs<T>(m_pRhsFun) );
    }

    virtual void localToGlobal_nonper(gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(gsMatrix<T>& globalRhs);

public: // *** Member functions ***

    virtual void assemble();

};

// ===================================================================================================================

/// @brief Visitor for right-hand side of the continuity equation (source function).
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSVisitorRhsP : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_localMat;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_trialFunData;
    using Base::m_periodicTransformMat;

public: // *** Constructor/destructor ***

    gsINSVisitorRhsP() {}
    
    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsINSVisitorRhsP(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr), m_pRhsFun(paramsPtr->getPde().source())
    {
        GISMO_ASSERT(m_pRhsFun == NULL || m_pRhsFun->targetDim() == 1, "Wrong RHS function passed into gsINSRhsP.");
    }
        

protected: // *** Member functions ***

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_trialUnkID = 1;    // pressure (not needed here)
    }

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_rhs<T>(m_pRhsFun) );
    }

    virtual void localToGlobal_nonper(gsMatrix<T>& globalRhs);

    virtual void localToGlobal_per(gsMatrix<T>& globalRhs)
    {
        localToGlobal_nonper(globalRhs);
    }

public: // *** Member functions ***

    virtual void assemble();

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSVisitors.hpp)
#endif