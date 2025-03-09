/** @file gsFlowTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gismo.h>


namespace gismo
{

template <class T>
class gsDistanceField
{
    
protected: // *** Class members ***

    gsMultiPatch<T> m_patches;
    gsMultiBasis<T> m_basis;
    gsBoundaryConditions<> m_bcDF;
    std::vector< gsField<T> > m_distanceField;
    index_t m_numRefs;
    std::vector<std::pair<int, boxSide> > m_bndIn;
    std::vector<std::pair<int, boxSide> > m_bndWall;


public: // *** Constructor/destructor ***

    gsDistanceField() {}
    
    gsDistanceField(gsMultiPatch<T> patches, gsMultiBasis<T> basis, index_t numRefs, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndWall) :
    m_patches(patches), m_basis(basis), m_numRefs(numRefs), m_bndIn(bndIn), m_bndWall(bndWall)
    {  }
    
    ~gsDistanceField() {}


public: // *** Class function ***

    gsField<T> computeDistanceField()
    {
        gsInfo << "Computing required distance field ... ";

        for (int i = 0; i < m_numRefs; ++i)       // additional refinements for distance computation
            m_basis.uniformRefine();

        gsFunctionExpr<real_t> fw("1", "0", "0", m_patches.targetDim());
        gsFunctionExpr<real_t> gw("0", "0", "0", m_patches.targetDim());
        gsFunctionExpr<real_t> wallw("0.0", "0.0", "0.0", m_patches.targetDim());

        // boundary conditions for related Poisson problem
        addBCs(m_bcDF, m_bndIn, m_bndWall, gw, wallw, 0);

        // solving Poisson problem
        gsPoissonAssembler<T> assembler(m_patches, m_basis, m_bcDF, fw, dirichlet::elimination, iFace::glue);
        assembler.assemble();

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

        // evaluating distance values at a grid of points
        size_t np = m_patches.nPatches();
        gsMultiPatch<T>* wallDistanceMP = new gsMultiPatch<T>;
        for (size_t i = 0; i < np; i++)
        {
            const gsBasis<T> & basis = m_basis.piece(i);

            std::vector< gsVector<T> > rr;
            rr.reserve(m_patches.parDim());

            for (short_t j = 0; j < m_patches.parDim(); ++j)            // computing grid of point
            {
                rr.push_back(basis.component(j).anchors().transpose());
            }
            gsMatrix<T> gridPts = gsPointGrid<T>(rr);

            unsigned geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            gsMapData<T> mapData;
            mapData.flags = geoFlags;
            mapData.patchId = i;
            mapData.points = gridPts;
            //typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *m_patches.patch(i)));
        
            /*
            gsMatrix<index_t> actives;
            gsMatrix<T> parGrads, physGrad;
            basis.deriv_into(gridPts, parGrads);
            //geoEval->evaluateAt(gridPts);

            for (index_t k = 0; k < gridPts.cols(); k++)
            {
                // eval solPoissonGrads at all pts
                basis.active_into(gridPts.col(k), actives);
                int numActP = actives.rows();
                gsMatrix<T> solActPoissonCoeffs(1, numActP);
                for (int j = 0; j < numActP; j++)
                    solActPoissonCoeffs.col(j) = solPoisson.coefficientVector(i).row(actives(j)).transpose();

                transformGradients(mapData, k, parGrads, physGrad);
                solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGrad.transpose();
            }
            */

            gsMatrix<T> solPoissonVals = solPoisson.value(gridPts, i);
            gsMatrix<T> solPoissonGrads = solPoisson.function(mapData.patchId).deriv(mapData.points);

            // Evaluate wall distance at pts
            gsMatrix<T> wallDistanceVals(1, gridPts.cols());
            for (index_t k = 0; k < gridPts.cols(); k++)
            {
                wallDistanceVals(0, k) = -solPoissonGrads.col(k).norm() + math::sqrt(math::pow(solPoissonGrads.col(k).norm(), 2) + 2 * solPoissonVals(0, k));
                if (math::isnan(wallDistanceVals(0, k)))
                {
                    wallDistanceVals(0, k) = 0.;
                }
            }
            
            typename gsGeometry<T>::uPtr geo = basis.interpolateAtAnchors(wallDistanceVals);    // interpolating distances at grid points 
            const gsMatrix<T> & distanceCoeffs = geo->coefs();
            wallDistanceMP->addPatch(basis.makeGeometry(distanceCoeffs));
        }

        gsInfo << "Done" << std::endl;

        return gsField<T>(m_patches, typename gsFunctionSet<T>::Ptr(wallDistanceMP), true);
    }
        

protected: // *** Class function ***

};


} // namespace gismo