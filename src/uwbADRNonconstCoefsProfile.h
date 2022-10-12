#pragma once

#include "uwbADRSolverLinear.h"
#include "uwbADRSolverNonlinear.h"
#include "uwbADRNonconstCoefsProfile.h"
#include "uwbRANSExamplesSetting.h"

namespace gismo
{

template<class T>
void solveADRprofile(const gsMultiPatch<T>& patches, std::vector<gsField<T> >& ADRcoeffs, uwbRANSProfileExample<real_t>& problemSettings)
{
    //---------------------------------------------------------------------------------------------
    /*gsVector<int> inInt(4);
    gsVector<real_t> inRealT(2);
    //gsVector<std::string> inString(3);
    gsVector<bool> inBool(7);

    std::ifstream inFile;
    inFile.open("ADRinitialSettingsProfile.txt");
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
    bool isoArtificialDiffusion = inBool(2);
    bool artificialDiffusion = inBool(3);
    bool FCT_lowOrder = inBool(4);
    bool animate = inBool(5);
    bool loadIC = inBool(6);

    int numElevate = inInt(0);
    int tauStabType = inInt(1);
    int numTimeIter = inInt(2);
    int numPicardIter = inInt(3);

    real_t delta_t = inRealT(0);
    real_t tol = inRealT(1);*/

    //----------------------------------------------------------------------------------------------

    bool supg = true;
    bool crosswind = true;
    bool isoArtificialDiffusion = false;
    bool artificialDiffusion = false;
    bool FCT_lowOrder = false;
    bool animate = false;
    bool loadIC = false;
    int numElevate = 0;
    int tauStabType = 1;
    int numTimeIter = 100;
    int numPicardIter = 5;
    real_t delta_t = 0.05;
    real_t tol = 0.001;

    //----------------------------------------------------------------------------------------------
    int plotPts = 10000;

    //------------- read from file ---------------------
    gsMatrix<real_t> ADRSolution;
    if (loadIC)
    {
        //gsFileData<> fdRead("ADR_ICsol.xml");
        gsFileData<> fdRead("ADRinitialConditionProfile.xml");
        ADRSolution = *(fdRead.getFirst< gsMatrix<real_t> >());
    }

    gsBoundaryConditions<> bcInfo;
    problemSettings.defineBCsADR(bcInfo, 1., 0.);

    // --------------- define Pde --------------

    uwbADRPde<real_t> ADRpde(patches, bcInfo);

    // --------------- set up basis ---------------
    // Copy basis from the geometry
    gsMultiBasis<> bases( patches );

    gsInfo << "degreeElevate = " << numElevate << "\n";
    bases.degreeElevate(numElevate);

    gsInfo << "bases.degree() = " << bases.degree() << "\n";

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
    params.settings().set(constantsADR::FCT_lowOrder, FCT_lowOrder); // use FCT low order stabilization
    params.settings().set(constantsADR::isoArtificialDiffusion, isoArtificialDiffusion); // use isotropic artificial diffusion stabilization
    params.settings().set(constantsADR::artificialDiffusion, artificialDiffusion); // use artificial diffusion stabilization
    if (supg || crosswind || artificialDiffusion)
        params.settings().set(constantsADR::tauStabType, tauStabType); //set formula for stabilization parameter tau

    params.settings().set(constantsADR::unst_innerIt, numPicardIter);

    //uwbADRSolverLinear<real_t> solverADR(params, diffusionCoeff, advectionCoeff, reactionCoeff);
    //necessary to use if crosswind or isotropic artificial diffusion stabilization is applied, even for linear ADR problem
    uwbADRSolverNonlinear<real_t> solverADR(params, ADRcoeffs[0], ADRcoeffs[1]);

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
    if(FCT_lowOrder)
        strStab += "FCTlowOrder_";
    if(supg)
        strStab += "SUPG_";
    if(crosswind)
        strStab += "CROSSWIND_";
    if(isoArtificialDiffusion)
        strStab += "ISOTROP_";
    if(artificialDiffusion)
        strStab += "ARTDIFF_";
    if (supg || crosswind || artificialDiffusion)
        strStab += "tau" + std::to_string(tauStabType);

    gsField<> velocity = solverADR.constructSolution();
    gsWriteParaview<>(velocity, "ADR_profile_" + strStab + "elev" + std::to_string(numElevate) + "Picard" + std::to_string(numPicardIter)
                      + "timeStep" + std::to_string(solverADR.getTimeStepNumber())
                      + "delta_t" + std::to_string(delta_t) + "tol" + std::to_string(tol), plotPts);

    gsWriteParaview<>(velocity, "ADRgrid_profile", plotPts, true);

    //--- save solution into file ---
    gsFileData<> fd;
    fd << solverADR.getSolution();
    fd.save("ADRsol_profile_" + strStab + "elev" + std::to_string(numElevate) + "Picard" + std::to_string(numPicardIter)
            + "timeStep" + std::to_string(solverADR.getTimeStepNumber())
            + "delta_t" + std::to_string(delta_t) + "tol" + std::to_string(tol) + ".xml");

}

} //namespace gismo
