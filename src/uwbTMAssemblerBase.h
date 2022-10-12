/** @file uwbTMAssemblerBase.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#ifdef _OPENMP 
#include <omp.h>
#endif

#include "uwbTMBlockAssembler.h"
#include "uwbTMEvaluators.h"

namespace gismo
{

template<class T>
class uwbINSSolverParams;

template<class T>
class uwbTMAssemblerBase
{
public:
    uwbTMAssemblerBase(uwbINSSolverParams<T>& params, int numVar) :
        m_blockAssembler(params, numVar)
    {
        std::string tmEvaluator = m_blockAssembler.getTMEvaluator();
        if (tmEvaluator == "koWilcoxLRN")
                m_pEvaluator = new uwbTMEvaluatorKOmegaWilcoxLRN<T>();
        else //(tmEvaluator == "koSST") or its variants
                m_pEvaluator = new uwbTMEvaluatorKOmegaSST<T>();

        initMembers();
    }

    virtual ~uwbTMAssemblerBase()
    {
    }

protected:
    void initMembers()
    {
        int numVar = getNumVar();
        int varDofs = numVarDofs();
        int dofs = numDofs();

        m_baseMatrix.resize(varDofs, dofs);

        m_matrix.resize(varDofs, dofs);
        m_rhs.setZero(varDofs, numVar);

        m_solution.setZero(varDofs, numVar);

        m_bInitialized = false;
        m_bSystemReady = false;
    }

public:
    virtual void initialize(const gsField<T> & uSolField)
    {
        initAssembly(uSolField);

        fillBase();

        m_bInitialized = true;
    }

    void setInitialCondition(const gsMatrix<T> & solVector)
    {
        m_solution = solVector;
        m_blockAssembler.setSolution(solVector);
    }

    virtual void update(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector, true);

        updateAssembly();

        fillSystem();
    }

    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_blockAssembler.updateCurrentSolField(solVector, false);

        updatePicardAssembly();

        fillSystem();
    }

protected:
    virtual void initAssembly(const gsField<T> & uSolField)
    { GISMO_NO_IMPLEMENTATION }

    virtual void updateAssembly()
    { GISMO_NO_IMPLEMENTATION }

    virtual void updatePicardAssembly()
    { GISMO_NO_IMPLEMENTATION }

    void checkReady() const
    { GISMO_ASSERT(m_bSystemReady, "System not ready, update() must be called first"); }

public:

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

    virtual void updateTurbCoeffSolField(const gsMatrix<T> & solVector, const gsField<T> & uSolField)
    { GISMO_NO_IMPLEMENTATION }

public:
    bool isInitialized() { return m_bInitialized; }
    const uwbTMBlockAssembler<T>& getBlockAssembler() const { return m_blockAssembler; }
    uwbTMBlockAssembler<T>& getBlockAssembler() { return m_blockAssembler; }
    int getNumVar() const { return m_blockAssembler.getNumVar(); }
    int numVarDofs() const { return m_blockAssembler.numVarDofs(); }
    int numDofs() const { return m_blockAssembler.numDofs(); }
    T getViscosity() const { return m_blockAssembler.getViscosity(); }

    const gsMatrix<T>&  getSolution() const { return m_solution; }
    virtual gsMatrix<T> getSolutionK_full(const gsDofMapper& pMapper, const gsMatrix<T>& solVector) const
    { GISMO_NO_IMPLEMENTATION }

    gsField<T> constructSolution(const gsMatrix<T>& solVector) const
    { return m_blockAssembler.constructSolution(solVector); }

public:

    void constructTurbCoeffSol(gsMultiPatch<T>& result, const gsMatrix<T> solution, const gsField<T>& uSolField, std::string coeffType = "turbViscosity")
    {
        result.clear(); // result is cleared first

        gsMatrix<T> turbCoeffVals;

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            gsMatrix<T> anchorPoints = m_blockAssembler.getBases().at(1).piece(p).anchors();

            evalTurbCoeff_singlePatch(p, anchorPoints, solution, uSolField, turbCoeffVals, coeffType);

            result.addPatch(m_blockAssembler.getBases().at(1).piece(p).interpolateAtAnchors(turbCoeffVals));
        }
    }

    gsField<T> constructTurbCoeffSol(const gsMatrix<T> solution, const gsField<T>& uSolField, std::string coeffType = "turbViscosity")
    {
        gsMultiPatch<T> * result = new gsMultiPatch<T>;
        constructTurbCoeffSol(*result, solution, uSolField, coeffType);
        return gsField<T>(m_blockAssembler.getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
    }

    virtual void evalTurbCoeff_into(index_t patchIndex, gsMatrix<T>& points, gsMatrix<T> solution, std::vector<gsMatrix<T> >& solUGrads, gsVector<T>& turbViscosity, std::string coeffType = "turbViscosity")
    { GISMO_NO_IMPLEMENTATION }

    gsField<T> constructTurbulentViscositySol(const gsMatrix<T> solution, const gsField<T>& uSolField)
    {
        return constructTurbCoeffSol(solution, uSolField, "turbViscosity");
    }

    void evalTurbViscosity_into(index_t patchIndex, gsMatrix<T>& points, gsMatrix<T> solution, std::vector<gsMatrix<T> >& solUGrads, gsVector<T>& turbViscosity)
    {
        evalTurbCoeff_into(patchIndex, points, solution, solUGrads, turbViscosity, "turbViscosity");
    }

protected:
    void evalTurbCoeff_singlePatch(index_t patchIndex, gsMatrix<T>& pts, const gsMatrix<T> solution, const gsField<T>& uSol, gsMatrix<T>& turbVals, std::string coeffType = "turbViscosity")
    {
        m_uSol = uSol;
        gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
        gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

        const int tarDim = patch->targetDim();

        std::vector<gsMatrix<T> > uGrads(pts.cols());//(np.prod());
        unsigned evFlags = NEED_GRAD_TRANSFORM;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> actives;
        gsMatrix<T> parGrads, physGrad;
        bases.front().deriv_into(pts, parGrads);
        geoEval->evaluateAt(pts);

        for (int k = 0; k < pts.cols(); k++)
        {
            // eval uGrads at all pts
            bases.front().active_into(pts.col(k), actives);
            int numActU = actives.rows();
            m_actCoeffsU.setZero(tarDim, numActU);

            for (int j = 0; j < numActU; j++)
                m_actCoeffsU.col(j) = uSol.coefficientVector(patchIndex).row(actives(j)).transpose();

            geoEval->transformGradients(k, parGrads, physGrad);
            uGrads[k].noalias() = m_actCoeffsU * physGrad.transpose();
        }

        // Evaluate turbulent viscosity at pts
        gsVector<T> turbCoeffVals;
        evalTurbCoeff_into(patchIndex, pts, solution, uGrads, turbCoeffVals, coeffType);
        turbVals = turbCoeffVals.transpose();
    }

public:
    void plotTurbulentViscosity(std::string const & fn, gsMatrix<T> solution, gsField<T> uSol, unsigned npts)
    {
        gsParaviewCollection collection(fn);
        std::string fileName;

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            fileName = fn + util::to_string(p);

            gsMatrix<T> geoVals, turbViscVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            evalTurbCoeff_singlePatch(p, pts, solution, uSol, turbViscVals, "turbViscosity");

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            if (3 - parDim > 0)
            {
                np.conservativeResize(3);
                np.bottomRows(3 - parDim).setOnes();
            }
            if (3 - tarDim > 0)
            {
                geoVals.conservativeResize(3, geoVals.cols());
                geoVals.bottomRows(3 - tarDim).setZero();
            }

            gsWriteParaviewTPgrid(geoVals, turbViscVals, np.template cast<index_t>(), fileName);

            collection.addPart(fileName, ".vts");
        }
        collection.save();
    }

    virtual uwbTMEvaluator<T>* getTMEvaluator() { return m_pEvaluator; }

protected:

    uwbTMBlockAssembler<T> m_blockAssembler;
    uwbTMEvaluator<T>* m_pEvaluator;

    bool m_bInitialized;
    bool m_bSystemReady;

    // base system matrix, which does not change at all
    gsSparseMatrix<T> m_baseMatrix;
    gsMatrix<T> m_baseRhs;

    // complete system matrix, nBaseMatrix + the part changing with Picard iteration
    gsSparseMatrix<T> m_matrix;
    gsMatrix<T> m_rhs;

    gsMatrix<T> m_solution;
    gsField<T> m_uSol;

    gsMatrix<T> m_actCoeffsU;

}; // class uwbTMAssemblerBase

//------------------------------------------------------------------------------------------------

template<class T>
class uwbTMAssemblerBaseUnsteady : public uwbTMAssemblerBase<T>
{
public:
    typedef uwbTMAssemblerBase<T> Base;

public:
    uwbTMAssemblerBaseUnsteady(uwbINSSolverParams<T>& params, int numVar) :
        Base(params, numVar)
    {
        m_blockAssembler.setUnsteady();
        m_timeStepSize = params.settings().get(constantsINS::timeStep);
        m_theta = params.settings().get(constantsINS::theta);

        initMembers();
    }

    virtual ~uwbTMAssemblerBaseUnsteady()
    {
    }

protected:
    void initMembers()
    {
        int varDofs = numVarDofs();

        m_nBaseMatrix.resize(varDofs, numDofs());
        m_nBaseRhs.setZero(varDofs, getNumVar());

        m_nExplicitPartRhs.setZero(varDofs, getNumVar());

        //---------------------------------------------------
        //TO DO: complete for theta method (stabilizations) and steady case!!!
        m_bAssembleAllBlocks = true;
        if (m_bAssembleAllBlocks)
            m_theta = 1;
        //---------------------------------------------------
    }

public:
    virtual void initialize(const gsField<T> & uSolField)
    {
        Base::initialize(uSolField);

        fillNBase();
    }

    virtual void update(const gsMatrix<T> & solVector, const gsField<T> & uSolField)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector, true);

        updateAssembly(uSolField);

        if (m_theta != 1)
        {
            fillExplicitPartMatrix();

            for (int s = 0; s < getNumVar(); s++)
                m_nExplicitPartRhs.col(s) = m_blockAssembler.getRhsF().col(s) - m_nExplicitPartMatrix.middleCols(s * this->numVarDofs(), this->numVarDofs()) * m_solution.col(s);
        }

        fillNBase();
        this->fillSystem();
    }

    virtual void fillNBase() { GISMO_NO_IMPLEMENTATION }
    virtual void fillExplicitPartMatrix() { GISMO_NO_IMPLEMENTATION }
    
    void changeTimeStepSize(const T timeStepSize)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_timeStepSize = timeStepSize;

        m_baseMatrix.resize(numVarDofs(), numDofs());
        this->fillBase();

        m_bSystemReady = false;
    }

    void setOldSolution(gsMatrix<T> oldSol) { m_oldSolution = oldSol; }
    void setTMblocksAssembly(bool assembleAllBlocks = true)
    {
        m_bAssembleAllBlocks = assembleAllBlocks;
        GISMO_ASSERT(m_theta == 1., "\nAll blocks assembly not ready for theta-method. Theta = 1 is set.\n");
        m_theta = 1.;
    }

protected:

    virtual void updateAssembly(const gsField<T> & uSolField)
    { GISMO_NO_IMPLEMENTATION }

protected:

    T m_timeStepSize;
    T m_theta;

    bool m_bAssembleAllBlocks;

    // base matrix + the part dependent on u_n, which does not change during one time step
    gsSparseMatrix<T> m_nBaseMatrix;
    gsMatrix<T> m_nBaseRhs;
    gsSparseMatrix<T> m_nExplicitPartMatrix;
    gsMatrix<T> m_nExplicitPartRhs;

    gsMatrix<T> m_oldSolution;

    // members from uwbTMAssemblerBase
    using Base::m_bInitialized;
    using Base::m_baseMatrix;
    using Base::m_bSystemReady;
    using Base::m_blockAssembler;
    using Base::m_solution;

public:
    // functions from uwbTMAssemblerBase
    using Base::numDofs;
    using Base::numVarDofs;
    using Base::getNumVar;
    using Base::getBlockAssembler;

}; // class uwbTMAssemblerBaseUnsteady

} //namespace gismo
