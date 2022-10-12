/** @file uwbINSSolversExample.cpp

Author(s): H. Hornikova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

// solvers
#include "uwbINSSolverSteady.h"
#include "uwbINSSolverUnsteady.h"
#include "uwbINSSolverDecoupledInterface.h"

using namespace gismo;

template<class T> gsMultiPatch<T> BSplineStep2D(T const & a = 1, T const & b = 2, T const & a_in = 1);

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 

    bool plot = false;
    bool outFile = false; // generate output file
    bool dg = false; // use DGFEM to connect patches
    bool supg = false;
    int numRefine = 3;
    int plot_pts = 10000;

    real_t viscosity = 0.1;
    char method = 'C'; // C - coupled, P - projection
    real_t timeStep = 0.1; // time step for unsteady computation
    real_t alpha_u = 0.6; // relaxation parameter for velocity in projection method
    real_t alpha_p = 1; // relaxation parameter for pressure in projection method
    int numThreads = 1; // number of threads for assembly

    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("file", "Write info into .txt file", outFile);
    cmd.addSwitch("dg", "Use DGFEM to connect patches", dg);
    cmd.addInt("r", "uniformRefine",
        "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("s", "plotSamples",
        "Number of sample points to use for plotting", plot_pts);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (numRefine < 0)
    {
        gsInfo << "Number of refinements must be non-negative, quitting.\n";
        return -1;
    }

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    if (dg)
        opt.intStrategy = iFace::dg;
    else
        opt.intStrategy = iFace::glue;

    gsInfo << "Solving the backward step example.\n";

    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> f("0", "0", 2); // external force
    gsFunctionExpr<> P("0", 2); // pressure boundary condition (decoupled case)
    gsFunctionExpr<> Uin("(-4*(y-1.5)^2 + 1)", "0", 2); // inlet velocity
    gsFunctionExpr<> Uwall("0", "0", 2); // wall velocity

    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uin, 0);

    // for decoupled solvers
    if (method == 'P')
    {
        bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, P, 1);
        bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, P, 1);
    }

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    real_t a = 8;
    real_t b = 2;
    real_t a_in = 1;

    patches = BSplineStep2D<real_t>(a, b, a_in);

    gsInfo << patches << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> tbasis(patches);
    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    // ========================================= Define solver ========================================= 

    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    // --------- other options ---------
    params.setNumThreads(numThreads); // set the number of threads for assembly

    params.settings().set(constantsINS::timeStep, timeStep); 
    params.settings().set(constantsINS::alpha_u, alpha_u);
    params.settings().set(constantsINS::alpha_p, alpha_p);
    //params.settings().set(constantsINS::unst_innerIt, 5); // max number of inner Picard iter. for unstready
    //params.settings().set(constantsINS::unst_innerTol, 1e-4); // stopping tolerance for inner Picard iter. in unstready
    params.settings().set(constantsINS::dec_innerIt, 5); // max inner iter. in iterative decoupled solver
    params.settings().set(constantsINS::dec_innerTol, 1e-4); //stopping tolerance for inner iter. in iterative decoupled solver
    //params.settings().set(constantsINS::omega, 10); // set angular velocity in case of rotation (assumed around x-axis)
    //params.settings().setProjVersion(decoupled::proj1); // two versions of projection method: proj1, proj2
    //params.settings().setDecoupledMethod(decoupled::iterative); // treatment off-diagonal blocks in case of periodic conditions (iterative of coupled)
    //params.settings().set(constantsINS::SUPG, supg); // use SUPG stabilization
    // ---------------------------------

    // solvers
    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver
    //uwbINSSolverUnsteady<real_t> navStokes(params); // unsteady coupled solver
    //uwbINSSolverDecoupledInterface<real_t> navStokes(params); // interface for decoupled solvers
    
    // decoupled solvers can be used as unsteady with time step dt = alpha_u / (1 - alpha_u)

    // ========================================= Solving ========================================= 
    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

    gsInfo << "initialization...\n";



    navStokes.initialize(); // steady solver


    navStokes.solve(10, 1e-5);

    real_t Tassembly = navStokes.getAssemblyTime();
    real_t Tsolve = navStokes.getSolverSetupTime();
    real_t Tsetupsolve = navStokes.getSolveTime();


    gsInfo << "Assembly time:" << Tassembly << "\n";
    gsInfo << "Solve time:" << Tsolve << "\n";
    gsInfo << "Solver setup time:" << Tsetupsolve << "\n";
    
    if (outFile)
    {
        time_t t = std::time(0);   // get time now
        struct tm * now = std::localtime(&t);
        std::stringstream filename;
        filename << "step_" << (now->tm_mon + 1) << "-" << now->tm_mday << "_"
            << now->tm_hour << "-" << now->tm_min << ".txt";

        std::ofstream ofile;
        ofile.open(filename.str().c_str());
        ofile << "Solving backward step example:\n";
        ofile << "DG = " << dg << "\n";
        ofile << "viscosity = " << viscosity << "\n";
        ofile << "method = " << method << "\n";
        ofile << "inlet velocity = " << Uin << "\n";
        //ofile << "angular velocity = " << omega << "\n";
        /*if (method == 'P')
            ofile << "alpha_u = " << alpha_u << "\n";*/
        ofile << "number of iterations = " << navStokes.getIterationNumber() << "\n";
        ofile << "last relative solution change = " << navStokes.solutionChangeRelNorm() << "\n";
        ofile << "Assembly time:" << Tassembly << "\n";
        ofile << "Solve time:" << Tsolve << "\n";
        ofile << "Solver setup time:" << Tsetupsolve << "\n";
        ofile.close();
    }

    // Optionally plot solution in paraview
    if (plot) 
    {
        gsField<> velocity = navStokes.constructSolution(0);
        //gsField<> relVelocity = navStokes.constructSolution(0, true); // in case of rotation
        gsField<> pressure = navStokes.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "step_velocity", plot_pts, true);
        gsWriteParaview<>(pressure, "step_pressure", plot_pts);
    }

    return 0; 
}

template<class T> gsMultiPatch<T> BSplineStep2D(T const & a, T const & b, T const & a_in)
{
    gsMultiPatch<T> mp;

    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, 0.0, a, b / 2));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, b / 2, a, b));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-a_in, b / 2, 0.0, b));

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::west, 2, boundary::east);
    mp.addAutoBoundaries();

    return mp;
}
