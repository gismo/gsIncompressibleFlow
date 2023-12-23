
#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolverSteady.h>
#include <gsIncompressibleFlow/src/gsINSSolverUnsteady.h>
#include <gsIncompressibleFlow/src/gsINSUtils.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>

using namespace gismo;


void solveOld(gsINSSolverBase<real_t>& NSsolver, int maxIt, real_t tol, bool plot, int plotPts, std::string id);
void solveNew(gsINSSolverBase<real_t>& NSsolver, gsINSSolverParams<real_t>& params, int maxIt, real_t tol, bool plot, int plotPts, std::string id, bool unst);


int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 

    bool steady = true;
    bool unsteady = false;

    int deg = 1;
    int numRefine = 3;
    int maxIt = 10;
    int picardIt = 5;
    real_t viscosity = 0.1;
    real_t timeStep = 0.1;
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;

    bool plot = false;
    int plotPts = 10000;
    int numThreads = 1; 

    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem in a backward facing step (BFS) domain.");

    cmd.addSwitch("steady", "Solve steady problem with direct linear solver", steady);
    cmd.addSwitch("unsteady", "Solve unsteady problem with direct linear solver", unsteady);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    cmd.addInt("d", "deg", "B-spline degree for geometry representation", deg);
    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addInt("t", "nthreads", "Number of threads for parallel assembly", numThreads);
    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsInfo << "Solving Navier-Stokes problem in a backward facing step (BFS) domain.\n";
    gsInfo << "viscosity = " << viscosity << "\n";

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    real_t a = 8;
    real_t b = 2;
    real_t a_in = 1;

    patches = BSplineStep2D<real_t>(deg, a, b, a_in);

    gsInfo << patches << "\n";


    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
    gsFunctionExpr<> f("0", "0", 2); // external force

    defineBCs_step(bcInfo, bndIn, bndOut, bndWall, 2); // bcInfo, bndIn, bndOut, bndWall are defined here


    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);
    
    refineBasis_step(basis, numRefine, 0, 0, 0, 0, 2, a, b);
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)


    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsINSSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setInt("numThreads",numThreads);

    if (steady)
    {
        gsINSSolverSteady<real_t> NSsolver(params);

        gsInfo << "\nSolving the steady problem with direct linear solver.\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        solveOld(NSsolver, maxIt, tol, plot, plotPts, "steadyOld");
        solveNew(NSsolver, params, maxIt, tol, plot, plotPts, "steadyNew", false);
    }

    if (unsteady)
    {
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("maxIt_picard", picardIt);
        params.options().setReal("tol_picard", picardTol);

        gsINSSolverUnsteady<real_t> NSsolver(params);

        gsInfo << "\nSolving the unsteady problem with direct linear solver.\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        solveOld(NSsolver, maxIt, tol, plot, plotPts, "unsteadyOld");
    }

    return 0; 
}


void solveOld(gsINSSolverBase<real_t>& NSsolver, int maxIt, real_t tol, bool plot, int plotPts, std::string id)
{
    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    // gsField<> uSol = NSsolver.constructSolution(0);
    // gsWriteParaview(uSol, "initVel_old");

    NSsolver.solve(maxIt, tol);

    //NSsolver.solveStokes(); // solve the Stokes problem
    //NSsolver.solveGeneralizedStokes(maxIt, tol); // solve the time-dependent Stokes problem (only in unsteady solvers)

    gsInfo << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";

    //gsInfo << "\nblockN =\n" << NSsolver.getAssembler()->getBlockAssembler().getBlockN() << "\n";
    // gsInfo << "\nmatrix norm =\n" << NSsolver.getAssembler()->matrix().norm() << "\n";
    // gsInfo << "\nrhs norm =\n" << NSsolver.getAssembler()->rhs().norm() << "\n";

    if (plot) 
    {
        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);
 
        gsInfo << "Plotting in Paraview...";
        gsWriteParaview<>(velocity, "BFS_" + id + "_velocity", plotPts, true);
        gsWriteParaview<>(pressure, "BFS_" + id + "_pressure", plotPts);
        gsInfo << " done.\n";
    }
}

void solveNew(gsINSSolverBase<real_t>& NSsolver, gsINSSolverParams<real_t>& params, int maxIt, real_t tol, bool plot, int plotPts, std::string id, bool unst)
{
    if (unst)
        params.options().setSwitch("unsteady", true);

    gsSparseMatrix<real_t, RowMajor> matBase, mat;
    gsSparseMatrix<real_t, RowMajor> blockA, blockN, blockBt;
    gsMatrix<real_t> rhsUbase, rhsUnonlin, rhsB, rhs, sol, newSol;

    int nnzPerRowU = 1;
    for (int i = 0; i < 2; i++)
        nnzPerRowU *= 2 * params.getBases()[0].maxDegree(i) + 1;

    int nnzPerRowP = 1;
    for (int i = 0; i < 2; i++)
        nnzPerRowP *= 2 * params.getBases()[1].maxDegree(i) + 1;

    int udofs = NSsolver.getAssembler()->getUdofs();
    int pdofs = NSsolver.getAssembler()->getPdofs();
    int pshift = NSsolver.getAssembler()->getPshift();
    int dofs = pshift + pdofs;

    blockA.resize(pshift, pshift);
    blockA.reserve(gsVector<int>::Constant(pshift, nnzPerRowU));

    blockN.resize(pshift, pshift);
    blockN.reserve(gsVector<int>::Constant(pshift, nnzPerRowU));

    blockBt.resize(pshift, pdofs);
    blockBt.reserve(gsVector<int>::Constant(pshift, nnzPerRowU));

    matBase.resize(dofs, dofs);
    matBase.reserve(gsVector<int>::Constant(dofs, nnzPerRowU + nnzPerRowP));

    mat.resize(dofs, dofs);
    mat.reserve(gsVector<int>::Constant(dofs, nnzPerRowU + nnzPerRowP));

    rhsUbase.setZero(pshift, 1);
    rhsUnonlin.setZero(pshift, 1);
    rhsB.setZero(dofs, 1);
    rhs.setZero(dofs, 1);
    sol.setZero(dofs, 1);
    newSol.setZero(dofs, 1);

    std::vector<gsMatrix<> > ddofs = NSsolver.getAssembler()->getBlockAssembler().getDirichletDofs();

    gsINSVisitorUUlin<real_t> visitorA(params);
    gsINSVisitorUUnonlin<real_t> visitorN(params);
    gsINSVisitorPU<real_t> visitorBt(params);

    visitorA.initialize();
    visitorN.initialize();
    visitorBt.initialize();

    gsField<> uSol = NSsolver.getAssembler()->constructSolution(sol, 0);
    visitorN.setCurrentSolution(uSol);

    gsWriteParaview(uSol, "initVel_new");

    for(size_t p = 0; p < params.getPde().patches().nPatches(); p++)
    {
        const gsBasis<>* patchBasisU = &(params.getBases()[0].piece(p));

        index_t nBasesU = patchBasisU->size();

        visitorA.initOnPatch(p);
        visitorN.initOnPatch(p);
        visitorBt.initOnPatch(p);

        for(index_t i = 0; i < nBasesU; i++)
        {
            visitorA.evaluate(i);
            visitorA.assemble();
            visitorA.localToGlobal(ddofs, blockA, rhsUbase);

            visitorN.evaluate(i);
            visitorN.assemble();
            visitorN.localToGlobal(ddofs, blockN, rhsUnonlin);

            visitorBt.evaluate(i);
            visitorBt.assemble();
            visitorBt.localToGlobal(ddofs, blockBt, rhsB);
        }
    }

    blockA.makeCompressed();
    blockN.makeCompressed();
    blockBt.makeCompressed();

    for (index_t i = 0; i < blockA.rows(); ++i)
        for (typename gsSparseMatrix<real_t, RowMajor>::InnerIterator it(blockA, i); it; ++it)
            matBase.coeffRef(i, it.col()) += it.value();

    for (index_t i = 0; i < blockBt.rows(); ++i)
        for (typename gsSparseMatrix<real_t, RowMajor>::InnerIterator it(blockBt, i); it; ++it)
        {
            matBase.coeffRef(i, pshift + it.col()) += it.value();
            matBase.coeffRef(pshift + it.col(), i) += it.value();
        }
            
    // =====================================================

    gsSparseSolver<>::LU linSolver;

    if (unst)
    {
        
    }
    else
    {
        int iter = 0;
        real_t relNorm = std::numeric_limits<real_t>::infinity();

        while ((relNorm > tol) && (iter < maxIt))
        {
            gsInfo << "Iteration number " << iter + 1 << "...";

            mat = matBase;
            for (index_t i = 0; i < blockN.rows(); ++i)
                for (typename gsSparseMatrix<real_t, RowMajor>::InnerIterator it(blockN, i); it; ++it)
                    mat.coeffRef(i, it.col()) += it.value();

            rhs = rhsB;
            rhs.topRows(pshift) += rhsUbase + rhsUnonlin;

            linSolver.compute(mat);
            newSol = linSolver.solve(rhs);

            gsMatrix<> solChange = sol - newSol;
            relNorm = solChange.norm() / newSol.norm();

            gsInfo << " Solution change relative norm: " << relNorm << "\n";

            //gsInfo << "\nblockN =\n" << blockN.topLeftCorner(udofs, udofs) << "\n";
            // gsInfo << "\nmatrix norm =\n" << mat.norm() << "\n";
            // gsInfo << "\nrhs norm =\n" << rhs.norm() << "\n";

            sol = newSol;
            uSol = NSsolver.getAssembler()->constructSolution(sol, 0);
            visitorN.setCurrentSolution(uSol);

            blockN *= 0;
            rhsUnonlin.setZero();

            for(size_t p = 0; p < params.getPde().patches().nPatches(); p++)
            {
                const gsBasis<>* patchBasisU = &(params.getBases()[0].piece(p));

                index_t nBasesU = patchBasisU->size();

                visitorN.initOnPatch(p);

                for(index_t i = 0; i < nBasesU; i++)
                {
                    visitorN.evaluate(i);
                    visitorN.assemble();
                    visitorN.localToGlobal(ddofs, blockN, rhsUnonlin);
                }
            }

            iter++;
        }
    }

}