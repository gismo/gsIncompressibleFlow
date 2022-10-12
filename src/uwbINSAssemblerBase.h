/** @file uwbINSAssemblerBase.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbINSBlockAssembler.h"
#include "uwbINSSolverParams.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerBase
{
public:
    uwbINSAssemblerBase(uwbINSSolverParams<T>& params) :
        m_blockAssembler(params)
    {
    }

    virtual ~uwbINSAssemblerBase()
    {
    }

protected:
    void initMembers() { GISMO_NO_IMPLEMENTATION }

    virtual void reinitMembers() { initMembers(); }

public:
    virtual void initialize() { GISMO_NO_IMPLEMENTATION }

    virtual void update(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector, true);

        updateAssembly();

        fillMatrix();
        fillRhs();
    }

    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_blockAssembler.updateCurrentSolField(solVector, false);

        updatePicardAssembly();

        fillMatrix();
        fillRhs();
    }

    virtual void fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const
    {
        m_blockAssembler.fillStokesSystem_into(stokesMatrix, stokesRhs);
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk, bool relative = false) const
    {
        return m_blockAssembler.constructSolution(solVector, unk, relative);
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, bool relative = false) const
    {
        return m_blockAssembler.constructSolution(solVector, unk, patchNumbers, sides, relative);
    }

    virtual gsField<T> constructSolutionCombined(const gsMatrix<T>& solVector, gsVector<size_t> relPatches) const
    {
        return m_blockAssembler.constructSolutionCombined(solVector, relPatches);
    }

    virtual T computeFlowRate(int patch, boxSide side, gsMatrix<T> solution) const
    {
        return m_blockAssembler.computeFlowRate(patch, side, solution);
    }

    virtual T computeDimensionlessWallDistance(gsMatrix<T> solution, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        return m_blockAssembler.computeDimensionlessWallDistance(solution, distancePatches, distanceSides, viscosity, reynoldsNumber, uFreeStream, maxYplus, npts, print, estimate);
    }

    T computeAspectRatio(bool minAR = false) { return m_blockAssembler.computeAspectRatio(minAR); }

    void addPressureOutletCondition(int patch, boxSide side)
    {
        gsMatrix<unsigned> boundaryIndicesOutlet = m_blockAssembler.getBases().back().basis(patch).boundary(side);

        std::vector< gsMatrix< index_t > > boundaryDofsToEliminate;
        boundaryDofsToEliminate.resize(m_blockAssembler.getBases().back().nBases());
        for (size_t i = 0; i < boundaryDofsToEliminate.size(); ++i)
        {
            boundaryDofsToEliminate[i].setZero(0, 0);
            if ( i == static_cast<size_t>(patch) )
                boundaryDofsToEliminate[i] = boundaryIndicesOutlet;
        }

        m_blockAssembler.markDofsAsEliminatedZeros(boundaryDofsToEliminate, 1);

        reinitMembers();
    }

    void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    {
        m_blockAssembler.markDofsAsEliminatedZeros(boundaryDofs, unk);

        reinitMembers();
    }

    virtual const gsSparseMatrix<T> & matrix() const { GISMO_NO_IMPLEMENTATION }

    virtual const gsMatrix<T> & rhs() const { GISMO_NO_IMPLEMENTATION }

    virtual const gsSparseMatrix<T>& getVelocityMassMatrix()
    { 
        return m_blockAssembler.getBlockM();
    }

    virtual const gsSparseMatrix<T>& getPressureMassMatrix()
    {
        return m_blockAssembler.getBlockMp();
    }

    void preparePCDboundary(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType = 0)
    {
        m_blockAssembler.preparePCDboundary(bndIn, bndOut, bndWall, bcType);
    }

    virtual void fillPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        Ap = m_blockAssembler.getPressurePoissonMatrix(assembAp, lumping);

        gsSparseMatrix<T> ApforFp;
        if (assembAp == assembFp)
            ApforFp = Ap;
        else
            ApforFp = m_blockAssembler.getPressurePoissonMatrix(assembFp, lumping);

        if(Fp.nonZeros()) // the matrix is not empty
            Fp += getViscosity() * ApforFp + m_blockAssembler.getPressureConvectionMatrix();
        else 
            Fp = getViscosity() * ApforFp + m_blockAssembler.getPressureConvectionMatrix();

        m_blockAssembler.applyPCDboundaryConditions(Ap, Fp, bcType);
    }

    void getPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        m_blockAssembler.preparePCDboundary(bndIn, bndOut, bndWall, bcType);
        this->fillPCDblocks(Ap, Fp, bcType, assembAp, assembFp, lumping);
    }

    virtual void evalResiduum(const gsMatrix<T> & solVector, std::vector<T> & residuum)
    { GISMO_NO_IMPLEMENTATION }

protected:

    virtual void updateAssembly()
    {
        m_blockAssembler.assembleNonlinearPart();
    }

    virtual void updatePicardAssembly()
    { GISMO_NO_IMPLEMENTATION }

    virtual void fillBase() { GISMO_NO_IMPLEMENTATION }

    virtual void fillMatrix() { GISMO_NO_IMPLEMENTATION }

    virtual void fillRhs() { GISMO_NO_IMPLEMENTATION }

public:
    void plot2DVorticity(std::string const & fn, gsMatrix<T>& solution, unsigned npts)
    {
        gsParaviewCollection collection(fn);
        std::string fileName;

        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            fileName = fn + util::to_string(p);

            gsMatrix<T> geoVals, vorticityVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            evalVorticity_singlePatch(p, pts, uSolField, vorticityVals);

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

            gsWriteParaviewTPgrid(geoVals, vorticityVals, np.template cast<index_t>(), fileName);

            collection.addPart(fileName, ".vts");
        }
        collection.save();
    }

    void plotPressureCoefficient(std::string const & fn, gsMatrix<T>& solution,
                                   index_t referencePatchIndex, gsVector<T>& referencePoint, T freeStreamVelocity,
                                   T rho, unsigned npts)
    {
        gsParaviewCollection collection(fn);
        std::string fileName;

        //gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);
        gsField<T> pSolField = m_blockAssembler.constructSolution(solution, 1);

        gsMatrix<T> referenceSolPVal = pSolField.value(referencePoint, referencePatchIndex);

        for (unsigned int patchIndex = 0; patchIndex < m_blockAssembler.getPatches().nPatches(); ++patchIndex)
        {
            fileName = fn + util::to_string(patchIndex);

            gsMatrix<T> geoVals, pressCoeffVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            //----------
            //evalPressureCoeff_singlePatch(patchIndex, pts, uSolField, pSolField, pressCoeffVals);
            //gsMatrix<T> solUVals = uSolField.value(pts, patchIndex);
            gsMatrix<T> solPVals = pSolField.value(pts, patchIndex);

            index_t nEvalPoints = pts.cols();
            gsVector<T> pressureCoeffVals(nEvalPoints);
            // Evaluate pressure coefficient at pts
            for (int k = 0; k < nEvalPoints; k++)
            {
                T pCoeffK = (solPVals.coeff(k, 0) - referenceSolPVal(0, 0)) / rho * math::pow(freeStreamVelocity, 2);
                //T pCoeffK = (solPVals.coeff(k, 0) - referenceSolPVal(0, 0)) / rho * math::pow((solUVals.col(k)).norm(), 2);
                //T pCoeffK = solPVals.coeff(k, 0) / rho * math::pow((solUVals.col(k)).norm(), 2);
                pressureCoeffVals(k) = pCoeffK;
            }

            pressCoeffVals = pressureCoeffVals.transpose();

            //----------

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

            gsWriteParaviewTPgrid(geoVals, pressCoeffVals, np.template cast<index_t>(), fileName);

            collection.addPart(fileName, ".vts");
        }
        collection.save();
    }

protected:
    void evalVorticity_singlePatch(index_t patchIndex, gsMatrix<T>& pts, const gsField<T>& uSolField, gsMatrix<T>& vortVals)
    {
        gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
        gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

        const int tarDim = patch->targetDim();

        unsigned evFlags = NEED_GRAD_TRANSFORM;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> actives;
        gsMatrix<T> parGrads, physGrad;
        bases.front().deriv_into(pts, parGrads);
        geoEval->evaluateAt(pts);

        gsMatrix<T> actCoeffsU;
        gsMatrix<T> uGrads;
        index_t nQuPoints = pts.cols();
        gsVector<T> vorticityVals(nQuPoints);
        // Evaluate vorticity at pts
        for (int k = 0; k < nQuPoints; k++)
        {
            // eval uGrads at all pts
            bases.front().active_into(pts.col(k), actives);
            int numActU = actives.rows();
            actCoeffsU.setZero(tarDim, numActU);

            for (int j = 0; j < numActU; j++)
                actCoeffsU.col(j) = uSolField.coefficientVector(patchIndex).row(actives(j)).transpose();

            geoEval->transformGradients(k, parGrads, physGrad);
            uGrads.noalias() = actCoeffsU * physGrad.transpose();

            T vorK = uGrads.coeff(1, 0) - uGrads.coeff(0, 1);
            vorticityVals(k) = vorK;
        }

        vortVals = vorticityVals.transpose();
    }
    
public:
    void setSolution(gsMatrix<T> solVector)
    {
        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector, true);
    }

    bool isInitialized() { return m_bInitialized; }
    const uwbINSBlockAssembler<T>& getBlockAssembler() const { return m_blockAssembler; }
    uwbINSBlockAssembler<T>& getBlockAssembler() { return m_blockAssembler; }
    virtual const gsMatrix<T>& getSolution() const { return m_solution; }
    virtual int numDofs() const { return m_blockAssembler.numDofs(); }
    virtual int getUdofs() const { return m_blockAssembler.getUdofs(); }
    virtual int getPdofs() const { return m_blockAssembler.getPdofs(); }
    virtual int getPshift() const { return m_blockAssembler.getPshift(); }
    int getTarDim() const { return m_blockAssembler.getTarDim(); }
    real_t getViscosity() const { return m_blockAssembler.getViscosity(); }
    bool isRotation() const { return m_blockAssembler.isRotation(); }

protected:
    uwbINSBlockAssembler<T>     m_blockAssembler;
    gsMatrix<T>                 m_solution;
    bool                        m_bInitialized;

}; //uwbINSAssemblerBase

} //namespace gismo
