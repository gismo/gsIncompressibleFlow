/** @file uwbINSTestPruh.cpp

@brief Navier-Stokes example for uwb testing/experimental purposes.

Author(s): H. Hornikova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

#include "uwbINSSolverSteady.h"
#include "uwbINSSolverUnsteady.h"
#include "uwbINSSolverDecoupledInterface.h"
//#include "uwbRANSSolver.h"
//#include "uwbTMSolverKOmega.h"

using namespace gismo;

gsBoundaryConditions<> defineBCs(const gsFunction<> & Uin, gsFunctionExpr<> & Uwall, char method, gsMatrix<real_t> transformMatrix);
template<class T> void refineBasis_RANS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal);

int main(int argc, char *argv[])
{

    // ========================================= Settings ========================================= 

    bool plot = false;
    bool dg = false; // use DGFEM to connect patches
    int plot_pts = 50000;

    int numRefine = 1;
    int numRefineLocal = 2;

    char method = 'C'; // C - coupled, P - projection (decoupled solvers)
    decoupled::method decMethod = decoupled::iterative; // in case of method P
    decoupled::projection projVersion = decoupled::proj2; // in case of method P

    real_t viscosity = 0.01;
    real_t alpha_u = 0.5; // relaxation parameter for velocity in projection method
    real_t alpha_p = 1; // relaxation parameter for pressure in projection method
    real_t timeStep = 0.01; // time step for unsteady computation
    //int numThreads = 2; // number of threads for assembly

    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("dg", "Use DGFEM to connect patches", dg);
    cmd.addInt("r", "uniformRefine",
        "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("s", "plotSamples",
        "Number of sample points to use for plotting", plot_pts);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (numRefine<0)
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

    gsInfo << "Solving flow in the stationary area.\n";

    // ========================================= Define transform ========================================= 

    gsMatrix<real_t> transformMatrix(3, 3);

    const real_t phi = -(1. / 7)*EIGEN_PI;
    const real_t cos = math::cos(phi);
    const real_t sin = math::sin(phi);
    transformMatrix(0, 0) = 1;
    transformMatrix(0, 1) = 0;
    transformMatrix(0, 2) = 0;
    transformMatrix(1, 0) = 0;
    transformMatrix(1, 1) = cos;
    transformMatrix(1, 2) = sin;
    transformMatrix(2, 0) = 0;
    transformMatrix(2, 1) = -sin;
    transformMatrix(2, 2) = cos;

    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> f("0", "0", "0", 3); // external force
    gsFunctionExpr<> Uin(" 3*(-(100 / 9) * (sqrt(y^2+z^2)-0.7)^2 + 1)", "0", "0", 3);
    gsFunctionExpr<> Uwall("0", "0", "0", 3); 

    bcInfo = defineBCs(Uin, Uwall, method, transformMatrix);

    // ========================================= Define geometry ========================================= 

    std::string input(MOTOR_DATA_DIR "/uwb-pilsen/pruh_withoutTrailingEdge.xml");

    gsFileData<> fileData(input);

    gsMultiPatch<>::uPtr patches;
    if (fileData.has< gsMultiPatch<> >())
    {
        patches = fileData.getFirst< gsMultiPatch<> >();
    }
    else
    {
        gsWarn << "Input file doesn't have a gsMultiPatch inside.\n";
        return -1;
    }

    gsInfo << *patches << "\n";

    // ========================================= Define basis ========================================= 

    gsMultiBasis<> tbasis(*patches);

    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();

    //refineBasis_RANS(tbasis, numRefine, numRefineLocal);

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    // ========================================= Define solver ========================================= 

    uwbINSPde<real_t> NSpde(*patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    params.settings().set(constantsINS::timeStep, timeStep);
    params.settings().set(constantsINS::alpha_u, alpha_u);
    params.settings().set(constantsINS::alpha_p, alpha_p);
    params.settings().setDecoupledMethod(decMethod);
    params.settings().setProjVersion(projVersion);

    //params.settings().set(constantsINS::unst_innerIt, 5); // number of inner Picard iterations for unsteady solver

    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver
    //uwbINSSolverUnsteady<real_t> navStokes(params); // unsteady coupled solver
    //uwbINSSolverDecoupledInterface<real_t> navStokes(params); // interface for decoupled solvers

    //uwbRANSSolver<real_t> navStokes(params);

    // ========================================= Solving =========================================

    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";



    gsInfo << "initialization...\n";


    navStokes.initialize();




navStokes.solve(20, 1e-5);

real_t Tassembly = navStokes.getAssemblyTime();
real_t Tsolve = navStokes.getSolverSetupTime();
real_t Tsetupsolve = navStokes.getSolveTime();


gsInfo << "Assembly time:" << Tassembly << "\n";
gsInfo << "Solve time:" << Tsolve << "\n";
gsInfo << "Solver setup time:" << Tsetupsolve << "\n";

    // Optionally plot solution in paraview
    if (plot)
    {
        gsField<> velocity = navStokes.constructSolution(0);
        gsField<> pressure = navStokes.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "pruh_velocity", plot_pts, true);
        gsWriteParaview<>(pressure, "pruh_pressure", plot_pts);

        // plot solution on the blade
        /*gsVector<index_t> bPatches(2);
        bPatches << 1, 1;
        std::vector<boxSide> bSides;
        bSides.push_back(boundary::east);
        bSides.push_back(boundary::west);

        gsField<> bladeVelocity = navStokes.constructSolution(0, bPatches, bSides);
        gsField<> bladePressure = navStokes.constructSolution(1, bPatches, bSides);

        // Write solution to paraview files
        gsWriteParaview<>(bladeVelocity, "pruh_bVelocity", plot_pts);
        gsWriteParaview<>(bladePressure, "pruh_bPressure", plot_pts);*/

        //system("paraview pruh_velocity.pvd&");
    }

    return 0;
}

gsBoundaryConditions<> defineBCs(const gsFunction<> & Uin, gsFunctionExpr<> & Uwall, char method, gsMatrix<real_t> transformMatrix)
{
    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);

    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Uwall, 0);

    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
    
    bcInfo.addPeriodic(0, boundary::west, 0, boundary::east, 3);
    bcInfo.addPeriodic(2, boundary::west, 2, boundary::east, 3);

    // non-periodic
    /*bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);*/

    bcInfo.setTransformMatrix(transformMatrix);

    if (method == 'P')
    {
        gsFunctionExpr<real_t> P("0", 3);

        bcInfo.addCondition(2, boundary::back, condition_type::dirichlet, &P, 1);
    }

    return bcInfo;
}

template<class T> void refineBasis_RANS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal)
{
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    gsMatrix<> box(3, 2);

    std::vector<real_t> parArea;

    for (int j = 0; j < basis.dim(); j++)
        for (int i = 0; i < basis.nBases(); i++)
            parArea.push_back(basis.piece(i).component(j).getMaxCellLength());       

    for (int j = 0; j < numRefineLocal; j++)
    {
        for (int i = 0; i < parArea.size(); i++)
            gsInfo << parArea[i] << "\n";

        index_t k = 0;

        for (int i = 0; i < basis.nBases(); i++)
        {
            box << 0, parArea[k], 0, 0, 0, 0;
            basis.refine(i, box);
            box << 1 - parArea[k], 1, 0, 0, 0, 0;
            basis.refine(i, box);
            k++;
        }

        for (int i = 0; i < basis.nBases(); i++)
        {
            box << 0, 0, 0, parArea[k], 0, 0;
            basis.refine(i, box);
            box << 0, 0, 1 - parArea[k], 1, 0, 0;
            basis.refine(i, box);
            k++;
        }

        box << 0, 0, 0, 0, 1 - parArea[k], 1;
        basis.refine(0, box);
        k++;

        box << 0, 0, 0, 0, 0, parArea[k];
        basis.refine(1, box);
        box << 0, 0, 0, 0, 1 - parArea[k], 1;
        basis.refine(1, box);
        k++;

        box << 0, 0, 0, 0, 0, parArea[k];
        basis.refine(2, box);

        for (int i = 0; i < parArea.size(); i++)
            parArea[i] /= 2;
    }
}
