/** @file francisTurbineExample.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt);

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 

    // solvers
    bool steady = false;
    bool unsteady = false;

    // domain definition
    std::string inFile = FLOW_DATA_DIR "francisTurbine_domain.xml";
    int nBlades = 14; // given by the domain geometry in inFile

    // discretization settings
    int numRefine = 0;
    int wallRefine = 0;

    // problem parameters
    real_t viscosity = 0.01;
    real_t inVelMag = 1;
    real_t inVelAngle = 20; // in degrees
    real_t omega = 0; // angular velocity for rotation

    // solver settings
    int maxIt = 10;
    int picardIt = 5;
    int linIt = 50;
    real_t timeStep = 0.1;
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;
    real_t linTol = 1e-6;
    std::string precond = "MSIMPLER_FdiagEqual";
    bool stokesInit = false; // start unsteady problem from Stokes solution
    
    // output settings
    bool quiet = false;
    bool plot = false;
    int plotPts = 50000;
    bool plotMesh = false;
    std::string outPath = "";
    
    // ---------------------------------------------------------------------------------

    gsCmdLine cmd(".");

    cmd.addSwitch("steady", "Solve steady problem", steady);
    cmd.addSwitch("unsteady", "Solve unsteady problem", unsteady);
    
    cmd.addString("", "domain", "Full path to the input xml file containing the domain geometry", inFile);
    cmd.addInt("b", "nBlades", "Number of blades (for the given domain geometry)", nBlades);

    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "wallRefine", "Number of h-refinement steps near blades", wallRefine);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addReal("", "inVelMag", "Magnitude of inflow velocity", inVelMag);
    cmd.addReal("", "inVelAng", "Angle of inflow velocity (in degrees)", inVelAngle);
    cmd.addReal("o", "omega", "Angular velocity for rotating frame of reference", omega);

    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);
    cmd.addInt("", "linIt", "Max. number of GMRES iterations", linIt);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);
    cmd.addReal("", "linTol", "Tolerance for iterative linear solver", linTol);
    cmd.addString("p", "precond", "Preconditioner type (format: PREC_Fstrategy, PREC = {PCD, PCDmod, LSC, AL, SIMPLE, SIMPLER, MSIMPLER}, Fstrategy = {FdiagEqual, Fdiag, Fmod, Fwhole})", precond);
    cmd.addSwitch("stokesInit", "Set Stokes initial condition", stokesInit);

    cmd.addSwitch("quiet", "Supress (some) terminal output", quiet);
    cmd.addSwitch("plot", "Plot the final result in ParaView format", plot);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addSwitch("plotMesh", "Plot the computational mesh", plotMesh);
    cmd.addString("", "outPath", "Path to the output directory", outPath);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if ( !(steady || unsteady) )
        gsWarn << "All computation flags set to false - nothing will be computed.\nPlease select at least one of the flags: --steady, --unsteady\n\n";

    // ========================= Define problem (geometry, BCs, rhs, basis) =========================

    gsInfo << "Reading geometry from file:\n" << inFile << "\n\n";
    gsMultiPatch<> mp;
    gsReadFile<> mpFile(inFile, mp);

    // inflow velocity components
    real_t inVelAngleRad = - inVelAngle * (EIGEN_PI / 180); // in radians
    std::string cosIn = util::to_string(math::cos(inVelAngleRad));
    std::string sinIn = util::to_string(math::sin(inVelAngleRad));
    std::string magStr = "(" + util::to_string(inVelMag) + "/sqrt(x^2+y^2))*";
    std::string rotX = util::to_string(-omega) + "*y";
    std::string rotY = util::to_string(omega) + "*x";
    std::string xIn = rotX + "+" + magStr + "(-" + cosIn + "*x + " + sinIn + "*y)";
    std::string yIn = rotY + "+" + magStr + "(-" + sinIn + "*x - " + cosIn + "*y)";

    gsFunctionExpr<> Uin(xIn, yIn, "0", 3);
    //gsFunctionExpr<> Uwall("0", "0", "0", 3);
    gsFunctionExpr<> Ublade(rotX, rotY, "0", 3);
    gsFunctionExpr<> f("0", "0", "0", 3);

    // transformation matrix for periodic sides
    real_t phi = -(2. / nBlades)*EIGEN_PI;
    real_t cos = math::cos(phi);
    real_t sin = math::sin(phi);
    gsMatrix<real_t> transformMatrix(3, 3);
    transformMatrix(0, 0) = cos;
    transformMatrix(0, 1) = -sin;
    transformMatrix(0, 2) = 0;
    transformMatrix(1, 0) = sin;
    transformMatrix(1, 1) = cos;
    transformMatrix(1, 2) = 0;
    transformMatrix(2, 0) = 0;
    transformMatrix(2, 1) = 0;
    transformMatrix(2, 2) = 1;

    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, Ublade, 0);   
    bcInfo.addCondition(0, boundary::back, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::front, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::back, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(2, boundary::front, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(2, boundary::back, condition_type::dirichlet, Ublade, 0);
    bcInfo.addPeriodic(0, boundary::west, 0, boundary::east, 3);
    bcInfo.addPeriodic(2, boundary::west, 2, boundary::east, 3);
    bcInfo.setTransformMatrix(transformMatrix);
    // bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);
    // bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
    // bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uwall, 0);
    // bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);

    gsInfo << "Solving Navier-Stokes problem in the Francis turbine runner wheel.\n";
    gsInfo << mp;
    gsInfo << "viscosity = " << viscosity << "\n";
    gsInfo << "omega = " << omega << "\n";
    gsInfo << "source function = " << f << "\n";
    gsInfo << "blade velocity: " << Ublade << "\n";

    // ========================================= Define basis ========================================= 

    gsMultiBasis<> basis(mp); 

    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    if (wallRefine)
    {
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            refineFirstKnotSpan<3, real_t>(basis, wallRefine, p, 0);
            refineLastKnotSpan<3, real_t>(basis, wallRefine, p, 0);
        }
    }

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)

    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(mp, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setSwitch("quiet", quiet);
    params.options().setString("lin.solver", "iter");
    params.options().setInt("lin.maxIt", linIt);
    params.options().setReal("lin.tol", linTol);
    params.options().setString("lin.precType", precond);
    params.options().setReal("omega", omega);

    gsOptionList solveOpt;
    solveOpt.addInt("maxIt", "", maxIt);
    solveOpt.addInt("plotPts", "", plotPts);
    solveOpt.addReal("tol", "", tol);
    solveOpt.addSwitch("plot", "", plot);
    solveOpt.addSwitch("plotMesh", "", plotMesh);
    solveOpt.addSwitch("stokesInit", "", stokesInit);
    solveOpt.addString("id", "", "");

    if (steady)
    {
        solveOpt.setString("id", "steady");

        gsINSSolverSteady<real_t, ColMajor> solver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the steady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("lin.precType") << "\n";
        solveProblem(solver, solveOpt);
    }
    
    if (unsteady)
    {
        solveOpt.setString("id", "unsteady");
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("nonlin.maxIt", picardIt);
        params.options().setReal("nonlin.tol", picardTol);

        gsINSSolverUnsteady<real_t, ColMajor> solver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the unsteady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("lin.precType") << "\n";
        solveProblem(solver, solveOpt);
    }

    return 0;
}

template<class T, int MatOrder>
void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt)
{
    gsStopwatch clock;

    // ------------------------------------
    // prepare strings for output filenames

    std::string idStr = "francisTurbine_" + opt.getString("id");

    // ------------------------------------
    // solve problem

    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

    gsINSSolverUnsteady<T, MatOrder>* pSolver = dynamic_cast<gsINSSolverUnsteady<T, MatOrder>* >(&NSsolver);

    if (pSolver && opt.getSwitch("stokesInit"))
        pSolver->solveStokes();

    NSsolver.solve(opt.getInt("maxIt"), opt.getReal("tol"), 0); // the last argument = min. number of iterations

    real_t totalT = clock.stop();

    gsInfo << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";
    gsInfo << "Total solveProblem time:" << totalT << "\n\n";

    gsFlowLinSystSolver_iter<T, MatOrder, gsGMRes<T> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<T, MatOrder, gsGMRes<T> >* >( NSsolver.getLinSolver() );
    reportLinIterations(linSolverPtr);

    // ------------------------------------
    // plot

    if (opt.getSwitch("plot")) 
    {
        bool rot = NSsolver.getParams()->isRotation();
        int plotPts = opt.getInt("plotPts");

        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);

        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, idStr + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, idStr + "_pressure", plotPts);

        if (rot)
        {
            gsField<> relVelocity = NSsolver.constructSolution(0, true);
            gsWriteParaview<>(relVelocity, idStr + "_relVelocity", plotPts, opt.getSwitch("plotMesh"));
        }
    }
}