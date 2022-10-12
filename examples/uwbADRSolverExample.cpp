/** @file uwbADRSolverExample.cpp

Author(s): E. Turnerova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

// solvers
#include "uwbADRSolverLinear.h"
#include "uwbADRSolverNonlinear.h"

using namespace gismo;

int main(int argc, char *argv[])
{
    //---------------------------------------------------------------------------------------------
    /*gsVector<int> inInt(5);
    gsVector<real_t> inRealT(4);
    //gsVector<std::string> inString(3);
    gsVector<bool> inBool(4);

    std::ifstream inFile;
    inFile.open("ADRinitialSettingsLshapeStable.txt");
    if (!inFile) {
        gsInfo << "Unable to open file";
        exit(1); // terminate with error
    }

    for(int i = 0; i < inBool.rows(); i++)
        inFile >> inBool(i);
    for(int i = 0; i < inInt.rows(); i++)
        inFile >> inInt(i);
    for(int i = 0; i < inRealT.rows(); i++)
        inFile >> inRealT(i);

    // ========================================= Settings =========================================

    bool supg = inBool(0);
    bool crosswind = inBool(1);
    bool animate = inBool(2);
    bool loadIC = inBool(3);

    int numIncrease = inInt(0);
    int numElevate = inInt(1);
    int numInitUniformRefine = inInt(2);
    int tauStabType = inInt(3);
    int numTimeIter = inInt(4);

    real_t diffusionCoeff = inRealT(0);
    real_t reactionCoeff = inRealT(1);
    real_t delta_t = inRealT(2);
    real_t tol = inRealT(3);*/

    //----------------------------------------------------------------------------------------------

    bool supg = true;
    bool crosswind = true;
    bool animate = false;
    bool loadIC = false;
    int numIncrease = 0;
    int numElevate = 0;
    int numInitUniformRefine  = 4;
    int tauStabType = 1;
    int numTimeIter = 100;
    real_t diffusionCoeff = 0.0001;
    real_t reactionCoeff = 0.;
    real_t delta_t = 0.05;
    real_t tol = 0.001;

    //----------------------------------------------------------------------------------------------
    int plotPts = 10000;
    int numPicardIter = 5;

    gsVector<> advectionCoeff(2);
    advectionCoeff << 3/sqrt(13), -2/sqrt(13);
    //advectionCoeff << 0, 1;

    //------------- read from file ---------------------
    gsMatrix<real_t> ADRSolution;
    if (loadIC)
    {
        //gsFileData<> fdRead("ADR_ICsol.xml");
        gsFileData<> fdRead("ADRinitialConditionLshape.xml");
        ADRSolution = *(fdRead.getFirst< gsMatrix<real_t> >());
    }

    // --------------- read geometry from file ---------------

    //std::string fileSrc( "lshape2d_3patches_thb.xml" );
    //gsMultiPatch<real_t> patches;
    //gsReadFile<real_t>( fileSrc, patches);
    //patches.computeTopology();

    gsMultiPatch<> patches;

    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0, -1, 1, 0));
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(-1, -1, 0, 0));
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(-1, 0, 0, 1));

    patches.addInterface(0, boundary::west, 1, boundary::east);
    patches.addInterface(1, boundary::north, 2, boundary::south);
    patches.addAutoBoundaries();

    gsInfo << "The domain is a "<< patches <<"\n";

    // --------------- add bonudary conditions ---------------
    gsBoundaryConditions<> bcInfo;

    gsFunctionExpr<> Uin("1", 2);
    //gsFunctionExpr<> Uin("if( y<=-1, if( x <=-0.25, if(x >= -0.75, 1, 0), 0 ), 0 )", 2);
    //gsFunctionExpr<> Uin("if( x<=-1, if( y<=-0.25, if(y >= -0.75, 1, 0), 0 ), 0 )", 2);
    gsFunctionExpr<> Uwall("0", 2);

    /*bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);*/

    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);

    // --------------- define Pde --------------

    uwbADRPde<real_t> ADRpde(patches, bcInfo);

    // --------------- set up basis ---------------
    // Copy basis from the geometry
    gsMultiBasis<> bases( patches );

    gsInfo << "degreeIncrease = " << numIncrease << "\n";
    gsInfo << "degreeElevate = " << numElevate << "\n";
    bases.degreeIncrease(numIncrease);
    bases.degreeElevate(numElevate);

    gsInfo << "bases.degree() = " << bases.degree() << "\n";

    // Number of initial uniform refinement steps:
    for (int i = 0; i < numInitUniformRefine; ++i)
      bases.uniformRefine();

    /*//---
    gsMatrix<> box_v0(2, 2);
    box_v0 << 0, 0, 0.9, 1;
    bases.refine(0, box_v0);
    bases.refine(1, box_v0);
    gsMatrix<> box_u0(2, 2);
    box_u0 << 0.9, 1, 0, 0;
    gsMatrix<> box_u1(2, 2);
    box_u1 << 0, 0.1, 0, 0;
    bases.refine(0, box_u1);
    bases.refine(1, box_u0);
    bases.refine(2, box_u0);
    //---*/

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.intStrategy = iFace::glue;
    opt.dirValues = dirichlet::interpolation;
    //opt.dirValues = dirichlet::l2Projection;
    gsDofMapper mapper;
    bases.getMapper(opt.dirStrategy, opt.intStrategy, bcInfo, mapper, 0);
    // ------------------------------
    real_t sol = 0.;
    gsMatrix<> initialCondition(mapper.freeSize(), 1);
    initialCondition.col(0).setConstant(sol);


    uwbADRSolverParams<real_t> params(ADRpde, bases, opt);
    params.settings().set(constantsADR::timeStep, delta_t);

    params.settings().set(constantsADR::SUPG, supg); // use SUPG stabilization
    params.settings().set(constantsADR::CROSSWIND, crosswind); // use crosswind stabilization
    if (supg || crosswind)
        params.settings().set(constantsADR::tauStabType, tauStabType); //set formula for stabilization parameter tau

    params.settings().set(constantsADR::unst_innerIt, numPicardIter);

    //uwbADRSolverLinear<real_t> solverADR(params, diffusionCoeff, advectionCoeff, reactionCoeff);
    //necessary to use if crosswind stabilization is applied, even for linear ADR problem
    uwbADRSolverNonlinear<real_t> solverADR(params, diffusionCoeff, advectionCoeff, reactionCoeff);

    if (loadIC)
        solverADR.setInitialCondition(ADRSolution);
    else
        solverADR.setInitialCondition(initialCondition);


    // --------------- adaptive refinement loop ---------------

    gsInfo << "numDofs = " << solverADR.numDofs() << "\n\n";
    solverADR.initialize();
    if (animate)
        solverADR.solveWithAnimation(numTimeIter, 1, tol);
    else
        solverADR.solve(numTimeIter, tol);

    std::string strStab = "";
    if(supg)
        strStab += "SUPG_";
    if(crosswind)
        strStab += "CROSSWIND_";
    if (supg || crosswind)
        strStab += "tau" + std::to_string(tauStabType);

    gsField<> velocity = solverADR.constructSolution();
    gsWriteParaview<>(velocity, "ADR_Lshape_" + strStab + "refUn" + std::to_string(numInitUniformRefine)
                      + "incr" + std::to_string(numIncrease) + "elev" + std::to_string(numElevate)
                      + "diff" + std::to_string(diffusionCoeff) + "reaction" + std::to_string(reactionCoeff) + "timeStep" + std::to_string(solverADR.getTimeStepNumber())
                      + "delta_t" + std::to_string(delta_t) + "tol" + std::to_string(tol), plotPts);

    gsWriteParaview<>(velocity, "ADRgrid_Lshape", plotPts, true);

    //--- save solution into file ---
    gsFileData<> fd;
    fd << solverADR.getSolution();
    fd.save("ADRsol_Lshape_" + strStab + "refUn" + std::to_string(numInitUniformRefine)
            + "incr" + std::to_string(numIncrease) + "elev" + std::to_string(numElevate)
            + "diff" + std::to_string(diffusionCoeff) + "reaction" + std::to_string(reactionCoeff) + "timeStep" + std::to_string(solverADR.getTimeStepNumber())
            + "delta_t" + std::to_string(delta_t) + "tol" + std::to_string(tol)+".xml");

    return 0;
}
