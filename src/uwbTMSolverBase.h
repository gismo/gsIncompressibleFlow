/** @file uwbTMSolverBase.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbINSSolverParams.h"

namespace gismo
{

template<class T>
class uwbTMSolverBase
{

public:
    uwbTMSolverBase()
    {
        m_pAssembler = NULL;
        m_iterationNumber = 0;
        m_assembT = 0;
        m_solsetupT = 0;
        m_solveT = 0;
    }

    virtual ~uwbTMSolverBase()
    {
        if (m_pAssembler)
        {
            delete m_pAssembler;
            m_pAssembler = NULL;
        }
    }

protected:
    void pardisoSetup(typename gsSparseSolver<T>::PardisoLU& solver)
    {
        solver.setParam(7, 15);
        solver.setParam(9, 13);
        solver.setParam(12, 0);
    }

public:
    virtual void initialize(const gsField<T>& uSolField)
    {
        //if(isTMSUPG())
        //    getAssembler()->updateTurbCoeffSolField(m_solution, uSolField);

        m_clock.restart();
        getAssembler()->initialize(uSolField);
        m_assembT += m_clock.stop();

        updateVelocitySolution(uSolField);
        gsInfo << "TM initialized.\n";
    }

    void setInitialCondition(const gsMatrix<T> & solVector) 
    { 
        m_solution = solVector;
        getAssembler()->setInitialCondition(solVector);
        m_iterationNumber = 0;
    }

    void updateVelocitySolution(const gsField<T>& uSolField) { m_uSolField = uSolField; }

    virtual void nextIteration() { GISMO_NO_IMPLEMENTATION; }

    gsField<T> constructSolution() const
    {
        return getAssembler()->constructSolution(m_solution);
    }

    gsField<T> constructTurbulentViscositySol()
    {
        return getAssembler()->constructTurbulentViscositySol(m_solution, m_uSolField);
    }

    gsField<T> constructTurbulentViscositySol(gsField<T> uSolField)
    {
        return getAssembler()->constructTurbulentViscositySol(m_solution, uSolField);
    }

    gsField<T> constructCoefficientSol(std::string coeffType)
    {
        return getAssembler()->constructTurbCoeffSol(m_solution, m_uSolField, coeffType);
    }

    void plotTurbulentViscosity(std::string const & fn, unsigned npts = 10000)
    {
        gsInfo << "Plotting turbulent viscosity in Paraview...\n";
        getAssembler()->plotTurbulentViscosity(fn, m_solution, m_uSolField, npts);
    }

    T solutionChangeRelNorm(index_t var) const
    {
        GISMO_ASSERT(var < getAssembler()->getNumVar(), "Solution index higher than number of variables.");

        T relNorm;

        if (m_iterationNumber)
        {
            gsMatrix<T> solChangeVector = getAssembler()->getSolution().col(var) - m_solution.col(var);
            relNorm = solChangeVector.norm() / m_solution.col(var).norm();

        }
        else
        {
            relNorm = std::numeric_limits<T>::infinity();
        }

        return relNorm;
    }

    T solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsMatrix<T> solChangeVector = solOld - solNew;
        T relNorm = solChangeVector.norm() / solNew.norm();

        return relNorm;
    }

    virtual void dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsInfo << "     TM solution change relative norm: ";

        for (int i = 0; i < solOld.cols(); i++)
            gsInfo << solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";
        
        gsInfo << "\n";
    }


    void setPoissonSolution(const gsField<T>& solPoisson, const gsMultiPatch<T>& patches, const gsMultiBasis<T>& bases)
    {
        getAssembler()->getBlockAssembler().setPoissonSolution(solPoisson, patches, bases);
    }

    void setPoissonSolution(gsMultiPatch<T> patches, gsMultiBasis<T> bases, gsBoundaryConditions<T> bcInfo, gsFunctionExpr<T> f, bool plot = false, unsigned npts = 10000)
    {
        gsPoissonAssembler<T> assembler(patches, bases, bcInfo, f, dirichlet::elimination, iFace::glue);
        assembler.assemble();
        //gsInfo << "numDofsPoisson: " << assembler.numDofs() << "\n";

        #ifdef GISMO_WITH_PARDISO
        typename gsSparseSolver<T>::PardisoLU solver;
        pardisoSetup(solver);
        #else
        typename gsSparseSolver<T>::LU solver;
        #endif 

        gsMatrix<T> solVector = solver.compute( assembler.matrix() ).solve ( assembler.rhs() );

        for (int i = 0; i < solVector.rows(); i++)
        {
            if (solVector(i, 0) < 0)
                solVector(i, 0) = math::abs(solVector(i, 0));
        }

        gsField<T> solPoisson = assembler.constructSolution(solVector, 0);
        setPoissonSolution(solPoisson, patches, bases);

        if (plot)
            gsWriteParaview<T>(solPoisson, "PoissonSolution", npts);
    }

    void plotWallDistance(std::string const & fn = "wallDistance", unsigned npts = 10000)
    {
        plotWallDistance(getAssembler()->getBlockAssembler().getPoissonSolution(), fn, npts);
    }

    void plotWallDistance(const gsField<T>& solPoisson, std::string const & fn = "wallDistance", unsigned npts = 10000)
    {
        gsInfo << "Plotting wall distance from Poisson solution in Paraview...\n";

        const unsigned n = getAssembler()->getBlockAssembler().getPoissonPatches().nPatches();
        gsParaviewCollection collection(fn);
        std::string fileName;

        for (unsigned i = 0; i < n; ++i)
        {
            fileName = fn + util::to_string(i);

            gsMatrix<T> geoVals, wallDistanceVals;
            gsVector<unsigned> np;
            evalWD_singlePatch(i, npts, solPoisson, geoVals, wallDistanceVals, np);
            gsWriteParaviewTPgrid(geoVals, wallDistanceVals, np.template cast<index_t>(), fileName);

            collection.addPart(fileName, ".vts");
        }
        collection.save();
    }

protected:
    void evalWD_singlePatch(index_t patchIndex, unsigned npts, const gsField<T>& solPoisson, gsMatrix<T>& geoVals, gsMatrix<T>& wallDistanceVals, gsVector<unsigned>& np)
    {
        gsGeometry<T>* patch = &getAssembler()->getBlockAssembler().getPoissonPatches().patch(patchIndex);
        const gsMultiBasis<T> basisPoisson = getAssembler()->getBlockAssembler().getPoissonBasis();

        const int tarDim = patch->targetDim();
        const int parDim = patch->domainDim();

        gsMatrix<T> ab = patch->support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        np = uniformSampleCount(a, b, npts);      
        gsMatrix<T> pts = gsPointGrid(a, b, np);

        std::vector<gsMatrix<T> > solPoissonGrads(np.prod());
        unsigned evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> actives;
        gsMatrix<T> parGrads, physGrad;
        basisPoisson.basis(patchIndex).deriv_into(pts, parGrads);
        geoEval->evaluateAt(pts);

        for (int k = 0; k < pts.cols(); k++)
        {
            // eval solPoissonGrads at all pts
            basisPoisson.basis(patchIndex).active_into(pts.col(k), actives);
            int numActP = actives.rows();
            gsMatrix<T> solActPoissonCoeffs(1, numActP);

            for (int j = 0; j < numActP; j++)
                solActPoissonCoeffs.col(j) = solPoisson.coefficientVector(patchIndex).row(actives(j)).transpose();

            geoEval->transformGradients(k, parGrads, physGrad);
            solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGrad.transpose();
        }

        gsMatrix<T> solPoissonVals = solPoisson.value(pts, patchIndex);

        // Evaluate geometry at pts
        geoVals = patch->eval(pts);

        // Evaluate wall distance at pts
        gsVector<T> wdVals(pts.cols());
        for (int k = 0; k < pts.cols(); k++)
        {
             wdVals(k) = -solPoissonGrads[k].norm() + math::sqrt(math::pow(solPoissonGrads[k].norm(), 2) + 2 * solPoissonVals(k));
            if (math::isnan( wdVals(k)))
            {
               wdVals(k) = 0.;
            }
        }
        wallDistanceVals = wdVals.transpose();

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
    }

public:
    void evalTurbViscosity_into(index_t patchIndex, gsMatrix<T>& points, std::vector<gsMatrix<T> >& solUGrads, gsVector<T>& turbViscosity) //const
    { getAssembler()->evalTurbViscosity_into(patchIndex, points, m_solution, solUGrads, turbViscosity); }

    virtual uwbTMAssemblerBase<T>*  getAssembler() const { return m_pAssembler; }
    virtual gsMatrix<T>     getSolutionK_full(const gsDofMapper& pMapper) const { return getAssembler()->getSolutionK_full(pMapper, m_solution); }
    const gsMatrix<T> &     getSolution() const { return m_solution; }
    bool isInitialized() { return getAssembler()->isInitialized(); }
    bool isTMSUPG() { return getAssembler()->getBlockAssembler().isTMsupg(); }
    std::string getTMEvaluator() { return getAssembler()->getBlockAssembler().getTMEvaluator(); }
    const T getAssemblyTime() const { return m_assembT; }
    const T getSolverSetupTime() const { return m_solsetupT; }
    const T getSolveTime() const { return m_solveT; }

protected:
    uwbTMAssemblerBase<T>* m_pAssembler;
    gsMatrix<T>            m_solution;
    gsField<T>             m_uSolField;
    unsigned               m_iterationNumber;

    gsStopwatch m_clock;
    T m_assembT, m_solsetupT, m_solveT;

}; // class uwbTMSolverBase

template<class T>
class uwbTMSolverBaseUnsteady : public uwbTMSolverBase<T>
{
public:
    typedef uwbTMSolverBase<T> Base;

public:
    uwbTMSolverBaseUnsteady(uwbINSSolverParams<T>& params) : Base()
    {
        GISMO_ASSERT(!(params.getBCs().numPeriodic()), "Use interface for periodic boundary conditions in turbulence model!");

        m_innerIter = params.settings().get(constantsINS::turb_innerIt);
        m_innerFirstIter = params.settings().get(constantsINS::turb_innerFirstIt);
        m_innerTol = params.settings().get(constantsINS::turb_innerTol);
    }

    virtual ~uwbTMSolverBaseUnsteady()
    { }

public:
    void changeTimeStepSize(const T timeStepSize) { getAssembler()->changeTimeStepSize(timeStepSize); }

    virtual void nextTimeStep() { this->nextIteration(); }

    const T getSimulationTime() const { return (this->m_iterationNumber * getAssembler()->getTimeStepSize()); }

    virtual uwbTMAssemblerBaseUnsteady<T>*  getAssembler() const { return dynamic_cast<uwbTMAssemblerBaseUnsteady<T>*>(m_pAssembler); }

protected:
    int m_innerIter;
    int m_innerFirstIter;
    T m_innerTol;

    // members from uwbTMSolverBase
    using Base::m_pAssembler;
    using Base::m_solution;

}; // class uwbTMSolverBaseUnsteady

} //namespace gismo
