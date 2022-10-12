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
#include "uwbADRNonconstCoefs.h"

using namespace gismo;

template<class T> gsMultiPatch<T> BSplineLshape(int const & deg = 2);
template<class T> gsTensorBSpline<2, T> BSplineRect(int deg, const T llx = 0, const T lly = 0, const T a = 1, const T b = 1);

int main(int argc, char *argv[])
{
    //---------------------------------------------------------------------------------------------
    gsVector<int> inInt(7);
    gsVector<real_t> inRealT(4);
    //gsVector<std::string> inString(3);
    gsVector<bool> inBool(7);

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
    bool isoArtificialDiffusion = inBool(2);
    bool artificialDiffusion = inBool(3);
    bool FCT_lowOrder = inBool(4);
    bool animate = inBool(5);
    bool loadIC = inBool(6);

    int degGeometry = inInt(0);
    int numIncrease = inInt(1);
    int numElevate = inInt(2);
    int numInitUniformRefine = inInt(3);
    int tauStabType = inInt(4);
    int numTimeIter = inInt(5);
    int numPicardIter = inInt(6);

    real_t viscosity = inRealT(0);
    real_t reactionCoeff = inRealT(1);
    real_t delta_t = inRealT(2);
    real_t tol = inRealT(3);

    //----------------------------------------------------------------------------------------------

    /*bool supg = true;
    bool crosswind = true;
    bool isoArtificialDiffusion = false;
    bool artificialDiffusion = false;
    bool FCT_lowOrder = false;
    bool animate = false;
    bool loadIC = false;
    int degGeometry = 2;
    int numIncrease = 0;
    int numElevate = 0;
    int numInitUniformRefine  = 4;
    int tauStabType = 1;
    int numTimeIter = 100;
    int numPicardIter = 5;
    real_t diffusionCoeff = 0.0001;
    real_t reactionCoeff = 0.;
    real_t delta_t = 0.05;
    real_t tol = 0.001;*/

    //----------------------------------------------------------------------------------------------
    int plotPts = 10000;

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

    gsMultiPatch<> patches = BSplineLshape<real_t>(degGeometry);

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
    //bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
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
    params.settings().set(constantsADR::FCT_lowOrder, FCT_lowOrder); // use FCT low order stabilization
    params.settings().set(constantsADR::isoArtificialDiffusion, isoArtificialDiffusion); // use isotropic artificial diffusion stabilization
    params.settings().set(constantsADR::artificialDiffusion, artificialDiffusion); // use artificial diffusion stabilization
    if (supg || crosswind || artificialDiffusion)
        params.settings().set(constantsADR::tauStabType, tauStabType); //set formula for stabilization parameter tau

    params.settings().set(constantsADR::unst_innerIt, numPicardIter);


    std::vector<gsField<> > advectionDiffusionCoeffs = computeADCoeffs(patches, bases, viscosity);
    //uwbADRSolverLinear<real_t> solverADR(params, diffusionCoeff, advectionCoeff, reactionCoeff);
    //necessary to use if crosswind or isotropic artificial diffusion stabilization is applied, even for linear ADR problem
    uwbADRSolverNonlinear<real_t> solverADR(params, advectionDiffusionCoeffs[0], advectionDiffusionCoeffs[1]);

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
    gsWriteParaview<>(velocity, "ADR_Lshape_" + strStab + "refUn" + std::to_string(numInitUniformRefine)
                      + "degGeom" + std::to_string(degGeometry) + "incr" + std::to_string(numIncrease) + "elev" + std::to_string(numElevate)
                      + "Picard" + std::to_string(numPicardIter)
                      + "viscosity" + std::to_string(viscosity) + "timeStep" + std::to_string(solverADR.getTimeStepNumber())
                      + "delta_t" + std::to_string(delta_t) + "tol" + std::to_string(tol), plotPts);

    gsWriteParaview<>(velocity, "ADRgrid_Lshape", plotPts, true);

    //--- save solution into file ---
    gsFileData<> fd;
    fd << solverADR.getSolution();
    fd.save("ADRsol_Lshape_" + strStab + "refUn" + std::to_string(numInitUniformRefine)
            + "degGeom" + std::to_string(degGeometry) + "incr" + std::to_string(numIncrease) + "elev" + std::to_string(numElevate)
            + "Picard" + std::to_string(numPicardIter)
            + "viscosity" + std::to_string(viscosity) + "timeStep" + std::to_string(solverADR.getTimeStepNumber())
            + "delta_t" + std::to_string(delta_t) + "tol" + std::to_string(tol) + ".xml");

    return 0;
}

template<class T> gsMultiPatch<T> BSplineLshape(const int &deg)
{
    gsMultiPatch<T> mp;

    T a = 1;
    T b = 1;
    T c = 0;

    mp.addPatch(BSplineRect(deg, c, -b, a, b));
    mp.addPatch(BSplineRect(deg, -a, -b, a, b));
    mp.addPatch(BSplineRect(deg, -a, c, a, b));

    /*mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0., -b, a, 0.));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-a, -b, 0., 0.));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-a, 0., 0., b));*/

    mp.addInterface(0, boundary::west, 1, boundary::east);
    mp.addInterface(1, boundary::north, 2, boundary::south);

    //mp.computeTopology();
    mp.addAutoBoundaries();

    return mp;
}

template<class T>
gsTensorBSpline<2, T> BSplineRect(int deg, const T llx, const T lly, const T a, const T b) // llx - lower left x, lly - lower left y
{
    gsKnotVector<T> kv(0, 1, 0, deg + 1); // first, last, inter, mult_end

    int n = kv.size() - deg - 1;
    gsMatrix<T> coef(n*n, 2);

    switch (deg)
    {
    case 1:
    {
        coef << llx + 0, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + b,
            llx + a, lly + b;
        break;
    }
    case 2:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 2), lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 2),
            llx + (a / 2), lly + (b / 2),
            llx + a, lly + (b / 2),
            llx + 0, lly + b,
            llx + (a / 2), lly + b,
            llx + a, lly + b;
        break;
    }
    case 3:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 3), lly + 0,
            llx + (2. / 3) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 3),
            llx + (a / 3), lly + (b / 3),
            llx + (2. / 3) * a, lly + (b / 3),
            llx + a, lly + (b / 3),
            llx + 0, lly + (2. / 3) * b,
            llx + (a / 3), lly + (2. / 3) * b,
            llx + (2. / 3) * a, lly + (2. / 3) * b,
            llx + a, lly + (2. / 3) * b,
            llx + 0, lly + b,
            llx + (a / 3), lly + b,
            llx + (2. / 3) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    case 4:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 4), lly + 0,
            llx + (2. / 4) * a, lly + 0,
            llx + (3. / 4) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 4),
            llx + (a / 4), lly + (b / 4),
            llx + (2. / 4) * a, lly + (b / 4),
            llx + (3. / 4) * a, lly + (b / 4),
            llx + a, lly + (b / 4),
            llx + 0, lly + (2. / 4) * b,
            llx + (a / 4), lly + (2. / 4) * b,
            llx + (2. / 4) * a, lly + (2. / 4) * b,
            llx + (3. / 4) * a, lly + (2. / 4) * b,
            llx + a, lly + (2. / 4) * b,
            llx + 0, lly + (3. / 4) * b,
            llx + (a / 4), lly + (3. / 4) * b,
            llx + (2. / 4) * a, lly + (3. / 4) * b,
            llx + (3. / 4) * a, lly + (3. / 4) * b,
            llx + a, lly + (3. / 4) * b,
            llx + 0, lly + b,
            llx + (a / 4), lly + b,
            llx + (2. / 4) * a, lly + b,
            llx + (3. / 4) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    case 5:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 5), lly + 0,
            llx + (2. / 5) * a, lly + 0,
            llx + (3. / 5) * a, lly + 0,
            llx + (4. / 5) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 5),
            llx + (a / 5), lly + (b / 5),
            llx + (2. / 5) * a, lly + (b / 5),
            llx + (3. / 5) * a, lly + (b / 5),
            llx + (4. / 5) * a, lly + (b / 5),
            llx + a, lly + (b / 5),
            llx + 0, lly + (2. / 5) * b,
            llx + (a / 5), lly + (2. / 5) * b,
            llx + (2. / 5) * a, lly + (2. / 5) * b,
            llx + (3. / 5) * a, lly + (2. / 5) * b,
            llx + (4. / 5) * a, lly + (2. / 5) * b,
            llx + a, lly + (2. / 5) * b,
            llx + 0, lly + (3. / 5) * b,
            llx + (a / 5), lly + (3. / 5) * b,
            llx + (2. / 5) * a, lly + (3. / 5) * b,
            llx + (3. / 5) * a, lly + (3. / 5) * b,
            llx + (4. / 5) * a, lly + (3. / 5) * b,
            llx + a, lly + (3. / 5) * b,
            llx + 0, lly + (4. / 5) * b,
            llx + (a / 5), lly + (4. / 5) * b,
            llx + (2. / 5) * a, lly + (4. / 5) * b,
            llx + (3. / 5) * a, lly + (4. / 5) * b,
            llx + (4. / 5) * a, lly + (4. / 5) * b,
            llx + a, lly + (4. / 5) * b,
            llx + 0, lly + b,
            llx + (a / 5), lly + b,
            llx + (2. / 5) * a, lly + b,
            llx + (3. / 5) * a, lly + b,
            llx + (4. / 5) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    case 6:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 6), lly + 0,
            llx + (2. / 6) * a, lly + 0,
            llx + (3. / 6) * a, lly + 0,
            llx + (4. / 6) * a, lly + 0,
            llx + (5. / 6) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 6),
            llx + (a / 6), lly + (b / 6),
            llx + (2. / 6) * a, lly + (b / 6),
            llx + (3. / 6) * a, lly + (b / 6),
            llx + (4. / 6) * a, lly + (b / 6),
            llx + (5. / 6) * a, lly + (b / 6),
            llx + a, lly + (b / 6),
            llx + 0, lly + (2. / 6) * b,
            llx + (a / 6), lly + (2. / 6) * b,
            llx + (2. / 6) * a, lly + (2. / 6) * b,
            llx + (3. / 6) * a, lly + (2. / 6) * b,
            llx + (4. / 6) * a, lly + (2. / 6) * b,
            llx + (5. / 6) * a, lly + (2. / 6) * b,
            llx + a, lly + (2. / 6) * b,
            llx + 0, lly + (3. / 6) * b,
            llx + (a / 6), lly + (3. / 6) * b,
            llx + (2. / 6) * a, lly + (3. / 6) * b,
            llx + (3. / 6) * a, lly + (3. / 6) * b,
            llx + (4. / 6) * a, lly + (3. / 6) * b,
            llx + (5. / 6) * a, lly + (3. / 6) * b,
            llx + a, lly + (3. / 6) * b,
            llx + 0, lly + (4. / 6) * b,
            llx + (a / 6), lly + (4. / 6) * b,
            llx + (2. / 6) * a, lly + (4. / 6) * b,
            llx + (3. / 6) * a, lly + (4. / 6) * b,
            llx + (4. / 6) * a, lly + (4. / 6) * b,
            llx + (5. / 6) * a, lly + (4. / 6) * b,
            llx + a, lly + (4. / 6) * b,
            llx + 0, lly + (5. / 6) * b,
            llx + (a / 6), lly + (5. / 6) * b,
            llx + (2. / 6) * a, lly + (5. / 6) * b,
            llx + (3. / 6) * a, lly + (5. / 6) * b,
            llx + (4. / 6) * a, lly + (5. / 6) * b,
            llx + (5. / 6) * a, lly + (5. / 6) * b,
            llx + a, lly + (5. / 6) * b,
            llx + 0, lly + b,
            llx + (a / 6), lly + b,
            llx + (2. / 6) * a, lly + b,
            llx + (3. / 6) * a, lly + b,
            llx + (4. / 6) * a, lly + b,
            llx + (5. / 6) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    default:
        GISMO_ERROR("Degree not implemented.");
        break;
    }

    return gsTensorBSpline<2, T>(kv, kv, give(coef));
}
