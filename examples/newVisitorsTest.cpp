
#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverSteady.h>
#include <gsIncompressibleFlow/src/gsINSSolverUnsteady.h>
#include <gsIncompressibleFlow/src/gsINSUtils.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>

using namespace gismo;

template<class NSsolverType>
void solve(NSsolverType& NSsolver, int maxIt, real_t tol, bool plot, int plotPts, std::string id);

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 

    bool steady = false;
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
        gsINSSolverDirectSteady<real_t> NSsolver1(params);

        gsInfo << "\nSolving the steady problem with direct linear solver.\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        solve< gsINSSolverBase<real_t> >(NSsolver, maxIt, tol, plot, plotPts, "steadyOld");
        solve< gsINSSolver<real_t> >(NSsolver1, maxIt, tol, plot, plotPts, "steadyNew");
    }

    if (unsteady)
    {
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("maxIt_picard", picardIt);
        params.options().setReal("tol_picard", picardTol);

        gsINSSolverUnsteady<real_t> NSsolver(params);
        gsINSSolverDirectUnsteady<real_t> NSsolver1(params);

        gsInfo << "\nSolving the unsteady problem with direct linear solver.\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        solve< gsINSSolverBase<real_t> >(NSsolver, maxIt, tol, plot, plotPts, "unsteadyOld");
        solve< gsINSSolver<real_t> >(NSsolver1, maxIt, tol, plot, plotPts, "unsteadyNew");
    }

    return 0; 
}

//=============================================================================================

template<class NSsolverType>
void solve(NSsolverType& NSsolver, int maxIt, real_t tol, bool plot, int plotPts, std::string id)
{
    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    NSsolver.solve(maxIt, tol);

    gsInfo << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";

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
