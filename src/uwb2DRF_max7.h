/** @file uwbRANSrofile_steadyInitial_periodic.cpp

Author(s): H. Hornikova, E. Turnerova
*/
#pragma once

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

// solvers
#include "uwbINSSolverSteady.h"
#include "uwbRANSSolver.h"
#include "uwbTMSolverKOmegaLinSteady.h"
#include "uwbTMSolverKOmega.h"
#include "uwbTurbineUtils.h"
#include "uwbBladeProfile.h"
#include "uwbKaplanTurbineRunnerBlade.h"

//#include "../../motor/uwb-pilsen/uwbINSSolverUnsteadyIterative.h"
//#include "uwbRANSSolverIterative.h"

#include <sstream>

#include <gsOptimizer/gsQualityMeasure.h>
#include <math.h>
#include <string.h>
#include "../jku/gsMotorUtils.h"
#include "../../../src/gsModeling/gsCoonsPatch.h"

#include "uwb2DRFGeom_max.h"
#include "uwbObjectiveFunctionEvaluator.h"


using namespace gismo;

template <typename T>
class uwbOptFlowSolver
{

public:

uwbOptFlowSolver()
{
}

~uwbOptFlowSolver()
{
}

public:
    mutable T objF;
    real_t tolRelNorm, tolObjRelVal, tolSteady;
    T timeStep, weightLift, weightVelocity, weightPressure, weightEfficiency, viscosity, turbIntensity;
    mutable int RANSSize, KOSize, num_dofs, maxNumIterRANS, minNumIterRANS; 
    int refinement, minNumIter, numIterSteadyNS, numIterKOmegaSteady, numIterFirst, numIterAll, numIterGrad, nP, numOfOpens, numRefineLocal_u1_e, numKnot_u1_e, numRefineLocal_u2, numKnot_u2;
    int numRefine, numRefineUniformLocal_v, numUniformKnot_v, numRefineLocal_v, numKnot_v, numRefineLocal_u0, numKnot_u0, numRefineLocal_u1_s, numKnot_u1_s, numRefineLocalFirstKnot_v;
    bool velAVG, loftEfficiency;

    const real_t PI = 3.141592653589793238463;
    int num_blades = 4;
    int plot_pts = 30000;
    T angularVelocity = 2*PI*538.0/60.0;
    T gravity = 9.81;
    T etah = 0.945;
    T H = 9.88929;
    T Hn = H*etah;
    T density = 997.;
    int number_of_f_parts = 3;

//    precond::type precType = precond::LSC_Adiag; //int precType = 1;
//    precond::type precType = precond::MSIMPLER_AdiagEqual;
//    int linMaxIt = 100;
//    T linTol = 1e-7;
//    T gamma = 2.5;
    int maxPicardIt = 5; 
    int maxTMPicardFirstIt = 50;
    std::string tmEvaluator;

    mutable gsVector<gsMatrix<T>> allOpreviousSolutionBest, allOinitObjParts, allOpTarget;
    gsVector<gsVector<T>> allOvOutTarget, allOviscosityRatio, allOvelocityRelativeX, allOvelocityRelativeY;
    gsVector<T> rr, allOflowRate, allOdiffToOptOpen, allOuniformityParam;
    gsVector<std::string> allOnameOfOpen;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

template<typename Tem>
gsVector<Tem> flowSolver(int IOP, int bladeOpen, swarmMember<Tem> & fish, opening<Tem> & otevreni, gsBSpline<Tem> & dCSS, gsBSpline<Tem> & dCPS, Tem & objFValue, gsVector<Tem> & initialObjP, gsVector<Tem> & pTarg, bool first, bool plot, std::string outfile) const
{
    gsVector<Tem> u = fish.x;
    Tem camberX = u(IOP);
    Tem camberY = u(IOP+nP);
    Tem leadingAngle = u(IOP+2*nP);
    Tem trailingAngle = u(IOP+3*nP);
    Tem thicknessX = u(IOP+4*nP);
    Tem thicknessY = u(IOP+5*nP);
    Tem outputAngle = u(IOP+6*nP);
    Tem radius = u(IOP+7*nP);
    Tem angle = u(IOP+8*nP) + allOdiffToOptOpen(bladeOpen);
    Tem endingOffset = u(IOP+9*nP);
    Tem chordLength = u(IOP+10*nP);
    Tem rotationCenterX = u(IOP+11*nP);
    Tem rotationCenterY = u(IOP+12*nP);

    gsVector<Tem> initSolution, viscosityRatio (nP);
    Tem uTangentialTarget, uniformity_param, velocity_absolute_x, velocity_absolute_y;

    Tem omega = 0.0;//2*PI*538.0/60.0;   // value 538 cycles per minute is given
    Tem velocity_blade = -rr(IOP)*omega;

    uTangentialTarget = allOvOutTarget(bladeOpen)(IOP); 
    viscosityRatio = allOviscosityRatio(bladeOpen);
    velocity_absolute_x = allOvelocityRelativeX(bladeOpen)(IOP);
    velocity_absolute_y = allOvelocityRelativeY(bladeOpen)(IOP) + velocity_blade;
    uniformity_param = allOuniformityParam(bladeOpen); 

    gsAssemblerOptions optAs;
    optAs.dirStrategy = dirichlet::elimination;
    optAs.intStrategy = iFace::glue;

    gsVector<Tem> length_x1_vec(nP);
    length_x1_vec << -0.186559, -0.177974, -0.169389, -0.160645, -0.152061, -0.143476, -0.134891;
    length_x1_vec *= 2.;
    Tem length_x1 = length_x1_vec(IOP);
    Tem length_x2 = 0.433544;
 
    if (!first) { initSolution = otevreni.previousSolution.col(IOP); }

    gsBoundaryConditions<Tem> bcInfo; 
    defineBCs_NS(bcInfo, velocity_absolute_x, velocity_absolute_y, velocity_blade);
    gsFunctionExpr<Tem> f("0", "0", 2);

    gsMultiPatch<Tem> patches;
    patches = BSplineProfile2DBetweenPatch<Tem>(IOP, length_x1, length_x2, ((2*PI*rr(IOP))/num_blades), camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength, angle, rotationCenterX, rotationCenterY, uniformity_param);

    patches.addInterface(0, boundary::north, 0, boundary::south);
    patches.addInterface(2, boundary::north, 2, boundary::south);
    patches.addAutoBoundaries();

    gsMultiBasis<Tem> tbasis(patches); // basis for RANS equations

    for (int i = 0; i < numRefine; ++i){tbasis.uniformRefine();}
    refineBasisUniformZones(tbasis, numRefineUniformLocal_v, numUniformKnot_v);
    if (refinement == 2){
        gsMatrix<> box_u0(2, 2);
        box_u0 << 0, 1, 0, 0;
        tbasis.refine(1, box_u0);
    }
refineBasisZones(tbasis, numRefineLocal_v, math::pow(2,numRefineUniformLocal_v) * numKnot_v, numRefineLocal_u0, numRefineLocal_u1_s, numRefineLocal_u1_e, numRefineLocal_u2, numKnot_u0, numKnot_u1_s, numKnot_u1_e, numKnot_u2); //for geometry_between = true
        //one more local refinement in v from 1st knot (in already refined mesh)
    refineBasisZones(tbasis, numRefineLocalFirstKnot_v, 1, 0, 0, 0, 0, 0, 0, 0, 0); //for geometry_between = true 

    //========================================================================================
    std::vector< gsMultiBasis<Tem> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    //================================= wall distance estimation ==================================================
    gsVector<Tem> Uin(2); Uin << velocity_absolute_x, velocity_absolute_y;
    Tem uFreeStream = Uin.norm();
    Tem Re = uFreeStream * (((2 * PI*rr(IOP)) / num_blades) / 2.0) / viscosity;
    gsVector<int> distancePatches(2); distancePatches << 1, 1;
    std::vector<boxSide> distanceSides;
    distanceSides.push_back(boundary::north);
    distanceSides.push_back(boundary::south);
    int numSamplePts = 50;
    Tem visc = viscosity;
     
    //============== construct NS =====================================
    Tem viscositySteady = 0.01;
    uwbINSPde<Tem> NSpde(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<Tem> params(NSpde, discreteBases, optAs);
    uwbINSSolverSteady<Tem> navStokes(params); // steady coupled solver
    navStokes.initialize(); // steady solver
    Tem maxYplus = 2.5;

    Tem wallDistance = navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, visc, Re, uFreeStream, maxYplus, numSamplePts, false, false);
    gsInfo << "wallDistance( " << IOP << "): " << wallDistance << "\n";

    // ========================================= Define turbulence solver ========================================= 
    gsBoundaryConditions<Tem> bcInfoTurb;
    Tem kInConst = 1.5 * math::pow(uFreeStream * turbIntensity, 2); // (3/2)*(UI)^2
    Tem kBlade = 1.5 * math::pow(velocity_blade * turbIntensity, 2);;
    Tem oInConst = kInConst / (viscosity * viscosityRatio(IOP)); // need to satisfy nu_T / nu approximately at inlet
    Tem oBlade;
    if (omega > 0.0) {oBlade = kBlade / (viscosity * viscosityRatio(IOP));}
    else {Tem beta = 0.0708; oBlade = 6 * viscosity / (beta * math::pow(wallDistance, 2));} 
    defineBCs_TM(bcInfoTurb, kInConst, kBlade, oInConst, oBlade);
    gsDofMapper koMapper;
    discreteBases[1].getMapper(optAs.dirStrategy, optAs.intStrategy, bcInfoTurb, koMapper, 0);

    if (first){ 
        //======================================= Init NS a K-omega ================================================
        uwbINSPde<Tem> koPdeSteady(patches, bcInfoTurb, f, viscositySteady);
        uwbINSSolverParams<Tem> koParamsSteady(koPdeSteady, discreteBases, optAs);
        gsMatrix<Tem> koInitial(koMapper.freeSize(), 2);
        koInitial.col(0).setConstant(kBlade);
        koInitial.col(1).setConstant(oBlade);
        uwbTMSolverKOmegaLinSteady<Tem> turbSolver(koParamsSteady);
        turbSolver.setInitialCondition(koInitial);

        //========================================== Solving steady NS =======================================
        navStokes.solve(numIterSteadyNS, tolSteady);
        RANSSize = navStokes.getSolution().rows();
        num_dofs = navStokes.numDofs();
        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

        // ========================================= Solving steady k-omega ========================================= 
        gsField<Tem> velocity = navStokes.constructSolution(0);
        turbSolver.initialize(velocity);
        turbSolver.solve(numIterKOmegaSteady, tolSteady); // solution change norm tol = 10^(-5)
        KOSize = turbSolver.getSolution().rows();
        gsVector<Tem> solStac (RANSSize+2*KOSize); 
        solStac.middleRows(0,RANSSize) = navStokes.getSolution();
        solStac.middleRows(RANSSize,KOSize) = turbSolver.getSolution().col(0);
        solStac.middleRows(RANSSize+KOSize,KOSize) = turbSolver.getSolution().col(1);
        initSolution = solStac;
        gsInfo << "TMnumDofs = " << koMapper.freeSize() << "\n";
        if (IOP==0){allOpreviousSolutionBest(bladeOpen).resize(RANSSize+2*KOSize,nP);}
    }

    // ========================================= Init RANS ========================================= 
    uwbINSPde<Tem> RANSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<Tem> RANSparams(RANSpde, discreteBases, optAs);

    RANSparams.settings().set(constantsINS::timeStep, timeStep);
    RANSparams.settings().set(constantsINS::unst_innerIt, maxPicardIt);

    uwbINSPde<Tem> koPde(patches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<Tem> koParams(koPde, discreteBases, optAs);

    koParams.settings().set(constantsINS::timeStep, timeStep);
    koParams.settings().set(constantsINS::turb_innerFirstIt, maxTMPicardFirstIt);
    koParams.settings().setTurbulenceEvaluator(tmEvaluator);

    uwbTMSolverKOmega<Tem> turbSolver_unsteady(koParams);
    if (tmEvaluator == "koSST"){solvePoissonEquation(patches, turbSolver_unsteady, "", true, plot_pts);}
    RANSparams.settings().setTurbulenceSolver(&turbSolver_unsteady);

/*//======================== iterative solver ============================================
    RANSparams.settings().setPrecondType(precType);
    RANSparams.settings().set(constantsINS::precAL_gamma, gamma); //only for msimpler (maybe)
    RANSparams.settings().set(constantsINS::iter_maxIt, linMaxIt);
    RANSparams.settings().set(constantsINS::iter_tol, linTol);
    uwbRANSSolverIterative<Tem, uwbGMResRight<Tem> > ransSolver(RANSparams);
*///=======================================================================================

    uwbRANSSolver<Tem> ransSolver(RANSparams);

    gsVector<Tem> normalStream(2);
    normalStream << 0.0, 1.0;
    normalStream = normalStream / normalStream.norm();

    //sides at which is computed the lift and pressure
    std::vector<boxSide> profileSides;
    profileSides.push_back(boundary::north);
    profileSides.push_back(boundary::south);

    //patches at which is computed the lift (must have the same length as the vector of sides!!)
    gsVector<int> profilePatchNumbers(2);
    profilePatchNumbers << 1, 1;

    //sides at which is computed the velocity part of the objective function
    std::vector<boxSide> outSides;
    outSides.push_back(boundary::east);

    //patches at which is computed the velocity part of the objective function (must have the same length as the vector of outSides!!)
    gsVector<int> outPatchNumbers(1);
    outPatchNumbers << 2;

    //======================================= Solving RANS =========================================
    gsVector<Tem> initRANS (RANSSize); initRANS = initSolution.middleRows(0, RANSSize);
    gsMatrix<Tem> initKO (KOSize, 2); 
    initKO.col(0) = initSolution.middleRows(RANSSize,KOSize);
    initKO.col(1) = initSolution.middleRows(RANSSize+KOSize,KOSize);

    turbSolver_unsteady.setInitialCondition(initKO);
    ransSolver.setInitialCondition(initRANS);
    ransSolver.initialize();

    gsField<Tem> uSol = ransSolver.constructSolution(0);
    gsField<Tem> pSol = ransSolver.constructSolution(1);

    uwbObjectiveFunctionEvaluator<Tem> objFunction(uSol, pSol, discreteBases, patches);
    objFunction.initialize(profilePatchNumbers, profileSides, normalStream, rr(IOP), outPatchNumbers, outSides, uTangentialTarget, velAVG);//false - compute velocity part by integral, true - by average
    objFunction.setWeights(weightLift, weightVelocity, weightPressure, weightEfficiency);  // weights for lift, velocity part at outlet, pressure distribution at profile, sum equal to 1!!!
    gsInfo << initialObjP << "\n";
 
    if (first) {objFunction.setTargetPressureAtProfile(-0.15,0.15);}// arbitrary values with no influence, objF is then nonsense
    else {objFunction.setTargetPressureAtProfile(pTarg); objFunction.setInitialObjectiveParts(initialObjP);}

    ransSolver.solveWithObjective(objFunction, tolObjRelVal, tolRelNorm, maxNumIterRANS, minNumIterRANS);
    
    otevreni.N_iter(IOP) = ransSolver.getIterationNumber();
    otevreni.CH(IOP) = ransSolver.getRelNorm();
    otevreni.CH_f(IOP) = ransSolver.getObjError();

    gsField<Tem> uSolution = ransSolver.constructSolution(0);
    gsField<Tem> pSolution = ransSolver.constructSolution(1);
    objFunction.updateSolution(uSolution, pSolution);

    if (first)
    {
        objFunction.setTargetPressureAtProfile(-0.15,0.15);
        pTarg = objFunction.getPressureTargetProfile();
    }
    gsVector<Tem> outputSolution (RANSSize+2*KOSize);
    outputSolution.middleRows(0,RANSSize) = ransSolver.getSolution();
    outputSolution.middleRows(RANSSize,KOSize) = turbSolver_unsteady.getSolution().col(0);
    outputSolution.middleRows(RANSSize+KOSize,KOSize) = turbSolver_unsteady.getSolution().col(1);
    objFValue = objFunction.evaluate();

    if (first){initialObjP = objFunction.getObjectives(); gsInfo << initialObjP << "\n";} 
    else {otevreni.f_parts.col(IOP) = objFunction.getObjectives();}

    gsInfo << "objFValue (blade " << IOP << "): " << objFValue << ", pTarg = " << pTarg << "\n";

    //================================================= Plot =====================================================================

    if (plot){
        gsInfo << "Plotting in Paraview...\n";
        gsField<Tem> velocity = ransSolver.constructSolution(0);
        gsField<Tem> pressure = ransSolver.constructSolution(1);
        gsField<Tem> kOmega = turbSolver_unsteady.constructSolution();

        gsWriteParaview<Tem>(velocity, outfile+"_velocity", plot_pts, false);
        gsWriteParaview<Tem>(pressure, outfile+"_pressure", plot_pts, false);
        gsWriteParaview<Tem>(kOmega, outfile+"_komega", plot_pts, true);
        ransSolver.plotTurbulentViscosity(outfile+"_turb_viscosity", plot_pts);
    }

    //======save pressure for 3D with loft=========
    const gsTensorBSplineBasis<2, T>*  basisBladePatch = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&(tbasis.basis(1)));
    gsBSpline<T> pSol_spline_pressureSide(basisBladePatch->knots(0),pSolution.coefficientVector(1).block(0,0,tbasis.at(1).component(0).size(),1));
    gsBSpline<T> pSol_spline_suctionSide(basisBladePatch->knots(0),pSolution.coefficientVector(1).block(tbasis.at(1).size()-tbasis.at(1).component(0).size()-1,0,tbasis.at(1).component(0).size(),1));

    dCSS = pSol_spline_suctionSide;
    dCPS = pSol_spline_pressureSide;

    return outputSolution;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
void initialize_vectors()
{
    allOpreviousSolutionBest.resize(numOfOpens);
    allOvOutTarget.resize(numOfOpens);
    allOinitObjParts.resize(numOfOpens);
    allOpTarget.resize(numOfOpens);
    allOflowRate.resize(numOfOpens);
    allOviscosityRatio = gsVector<gsVector<T>> (numOfOpens);//.resize(numOfOpens);
    allOvelocityRelativeX.resize(numOfOpens);
    allOvelocityRelativeY.resize(numOfOpens);
    allOdiffToOptOpen.resize(numOfOpens);
    allOuniformityParam.resize(numOfOpens);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
void flowSolver_load(int nPopulation, std::string text) const
{
    for (int nOpen = 0; nOpen < numOfOpens; nOpen++){
        allOpTarget(nOpen) = gsMatrix<T>::Zero(2,nP);
        allOinitObjParts(nOpen) = gsMatrix<T>::Zero (number_of_f_parts,nP);

        gsFileData<real_t> fdRead(allOnameOfOpen(nOpen) + "_" + text + "_RANS_solution.xml");
        gsMatrix<real_t> loadRANSSolution = *(fdRead.getFirst< gsMatrix<real_t> >());
        gsFileData<real_t> fdReadTM(allOnameOfOpen(nOpen) + "_" + text + "_TM_solution.xml");
        gsMatrix<real_t> loadTMSolution = *(fdReadTM.getFirst< gsMatrix<real_t> >());
        gsFileData<real_t> fdReadTP(allOnameOfOpen(nOpen) + "_" + text + "_targetPressure.xml");
        gsMatrix<real_t> TP = *(fdReadTP.getFirst< gsMatrix<real_t> >());
        gsFileData<real_t> fdReadTV(allOnameOfOpen(nOpen) + "_" + text + "_initObjParts.xml");
        gsMatrix<real_t> TV = *(fdReadTV.getFirst< gsMatrix<real_t> >());

        RANSSize = loadRANSSolution.rows();
        KOSize = loadTMSolution.rows();
                    
        for (int IOP = 0; IOP < nP; IOP++){
    
            gsVector<T> loadSolution (RANSSize + 2*KOSize);
            loadSolution.middleRows(0,RANSSize) = loadRANSSolution.col(IOP);
            loadSolution.middleRows(RANSSize,KOSize) = loadTMSolution.col(2*IOP);
            loadSolution.middleRows(RANSSize+KOSize,KOSize) = loadTMSolution.col(2*IOP+1);
            
            if (IOP == 0){
                allOpreviousSolutionBest(nOpen).resize(RANSSize + 2*KOSize,nP);
            }
    
            allOpreviousSolutionBest(nOpen).col(IOP) = loadSolution;
            allOpTarget(nOpen).col(IOP) = TP.col(IOP);
            allOinitObjParts(nOpen).col(IOP) = TV.col(IOP);    
        } 
    }
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void flowSolver_save(gsVector<T> u, std::string text, bool initValues) const
{
    gsVector<T> initSolution (RANSSize + 2*KOSize);
    gsMatrix<T> initRANS (RANSSize, nP), initKO (KOSize, 2*nP);
    gsMatrix<T> TP (2,nP), TV (number_of_f_parts,nP);
    gsMatrix<T> dP (u.size(),1); dP.col(0) = u;
    gsFileData<T> fd_dP; fd_dP << dP;
 
    fd_dP.save(text + "_profile_designParameters.xml");
 
    for (int nOpen = 0; nOpen < numOfOpens; nOpen++){
        for (int IOP = 0; IOP < nP; IOP++){
            initSolution = allOpreviousSolutionBest(nOpen).col(IOP); 
            TP.col(IOP) = allOpTarget(nOpen).col(IOP);
            TV.col(IOP) = allOinitObjParts(nOpen).col(IOP);

            initRANS.col(IOP) = initSolution.middleRows(0, RANSSize);
            initKO.col(2*IOP) = initSolution.middleRows(RANSSize,KOSize);
            initKO.col(2*IOP+1) = initSolution.middleRows(RANSSize+KOSize,KOSize);
        }    

        gsFileData<T> fd, fd_TM, fd_TP, fd_dP, fd_TV;
        fd << initRANS;
        fd_TM << initKO;    
        fd_TP << TP;
        fd_TV << TV;

        fd.save(allOnameOfOpen(nOpen) + "_" + text + "_RANS_solution.xml");
        fd_TM.save(allOnameOfOpen(nOpen) + "_" + text + "_TM_solution.xml"); 
        if (initValues){
            fd_TP.save(allOnameOfOpen(nOpen) + "_" + text + "_targetPressure.xml"); 
            fd_TV.save(allOnameOfOpen(nOpen) + "_" + text + "_initObjParts.xml");
        }
    }
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void flowSolver_draw(swarmMember<T> & fish, std::string text) const
{
    minNumIterRANS = 0;
    maxNumIterRANS = 0;
    for (int nOpen = 0; nOpen < numOfOpens; nOpen++){
        for (int IOP = 0; IOP < nP; IOP++){
            gsVector<T> init_f_parts = allOinitObjParts(nOpen).col(IOP);
            gsVector<T> p_targ = allOpTarget(nOpen).col(IOP);
            gsBSpline<T> dCSS, dCPS;
            T objfunct;
            std::ostringstream strs_profile;  strs_profile << IOP; 
            gsVector<T> empty = flowSolver<T>(IOP,nOpen,fish,fish.open[nOpen],dCSS,dCPS,objfunct,init_f_parts,p_targ,false,true,allOnameOfOpen(nOpen) + "_" + text + "_profile" + strs_profile.str());
        }
    }
    minNumIterRANS = minNumIter;
    maxNumIterRANS = numIterAll;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

T flowSolver_swarm(int nPopulation, swarmMember<T> & fish, bool start) const
{
    gsVector<T> objFsw_open (nP); objFsw_open.setZero(nP);
    T objFsw = 0.0;
    gsVector<std::vector<gsBSpline<T>>> defcurves_suctionside_pSol (numOfOpens), defcurves_pressureside_pSol (numOfOpens);    
    minNumIterRANS = minNumIter;

    for (int nOpen = 0; nOpen < numOfOpens; nOpen++){ 
        fish.open[nOpen].f = 0.0;
        std::vector<gsBSpline<T>> dCSSVec (nP), dCPSVec (nP); 

        if (start)
        {
            maxNumIterRANS = numIterFirst;
            allOinitObjParts(nOpen) = gsMatrix<T>::Zero(number_of_f_parts,nP); 
            allOpTarget(nOpen) = gsMatrix<T>::Zero(2,nP);                      

            #pragma omp parallel for num_threads(nP)
            for (int IOP = 0; IOP < nP; IOP++){
                gsVector<T> empty, init_f_parts, p_targ;
                gsBSpline<T> dCSS, dCPS;

                empty = flowSolver<T>(IOP, nOpen, fish, fish.open[nOpen], dCSS, dCPS, objFsw_open(IOP), init_f_parts, p_targ, true, false, "");
                allOpreviousSolutionBest(nOpen).col(IOP) = empty;
                allOinitObjParts(nOpen).col(IOP) = init_f_parts;
                allOpTarget(nOpen).col(IOP) = p_targ;
                dCSSVec[IOP] = dCSS;
                dCPSVec[IOP] = dCPS;
            }
        }
        else
        {
            maxNumIterRANS = numIterAll;
            //#pragma omp parallel for num_threads(7)
            for (int IOP = 0; IOP < nP; IOP++){
                gsVector<T> solR;
                gsVector<T> init_f_parts = allOinitObjParts(nOpen).col(IOP);
                gsVector<T> p_targ = allOpTarget(nOpen).col(IOP);
                gsBSpline<T> dCSS, dCPS;

                solR = flowSolver<T>(IOP, nOpen, fish, fish.open[nOpen], dCSS, dCPS, objFsw_open(IOP), init_f_parts, p_targ, false, false, "");
                fish.open[nOpen].previousSolution.col(IOP) = solR;
                dCSSVec[IOP] = dCSS;
                dCPSVec[IOP] = dCPS;
            }
        }

        defcurves_suctionside_pSol(nOpen) = dCSSVec;
        defcurves_pressureside_pSol(nOpen) = dCPSVec;

        for (int IOP = 0; IOP < nP; IOP++){
            fish.open[nOpen].f += objFsw_open(IOP)/nP;
        }
        objFsw += fish.open[nOpen].f;
    }

    if (loftEfficiency){

        //======pressure 3D=========
        gsKnotVector<T> kvloft(0, 1, nP-4, 4);
    
        gsTensorBSpline<2,T> suction_side_surface_pSol;
        gsTensorBSpline<2,T> pressure_side_surface_pSol;
        gsTensorBSpline<2,T> suction_side_surface_blade3D;
        gsTensorBSpline<2,T> pressure_side_surface_blade3D;
    
        gsMatrix<T> section_parameters(1,nP);
        for (index_t i = 0; i < nP; i++) {
            section_parameters(0,i) = (rr[i]-rr[0])/(rr[nP-1]-rr[0]);
        }

        gsVector<T> A, posun(nP);

        for (int nOpen = 0; nOpen < numOfOpens; nOpen++){
            computeLoftFuction(defcurves_suctionside_pSol(nOpen), defcurves_suctionside_pSol(nOpen)[0].basis().knots(), nP, kvloft, section_parameters, suction_side_surface_pSol, 1);
            computeLoftFuction(defcurves_pressureside_pSol(nOpen), defcurves_pressureside_pSol(nOpen)[0].basis().knots(), nP, kvloft, section_parameters, pressure_side_surface_pSol, 1);
    
            posun.setConstant(nP,allOdiffToOptOpen(nOpen));
            A = fish.x.middleRows(8*nP,nP) + posun;
            bool extrapolateBlade = false;
            KaplanTurbineRunnerBlade<T> runnerBlade(nP, fish.x.middleRows(0*nP,nP), fish.x.middleRows(1*nP,nP), fish.x.middleRows(2*nP,nP), fish.x.middleRows(3*nP,nP), fish.x.middleRows(4*nP,nP),          
                                                    fish.x.middleRows(5*nP,nP), fish.x.middleRows(9*nP,nP), fish.x.middleRows(6*nP,nP), fish.x.middleRows(7*nP,nP), fish.x.middleRows(10*nP,nP),
                                                    A, fish.x.middleRows(11*nP,nP), fish.x.middleRows(12*nP,nP), rr, extrapolateBlade);

            runnerBlade.compute(false, false);

            pressure_side_surface_blade3D = runnerBlade.getPressureSideSurface();
            suction_side_surface_blade3D = runnerBlade.getSuctionSideSurface();

            gsField<T> pSol_3D_suction(suction_side_surface_blade3D, suction_side_surface_pSol);
            gsField<T> pSol_3D_pressure(pressure_side_surface_blade3D, pressure_side_surface_pSol);

            const gsBasis<T>& basisSuction = suction_side_surface_pSol.basis();
            const gsBasis<T>& basisPressure = pressure_side_surface_pSol.basis();

            //const gsBasis<T>& basisInlet = surface_pIn.basis();
            //const gsBasis<T>& basisOutlet = surface_pOut.basis();

            int numBlades = num_blades;
            uwbObjectiveFunctionEvaluator<T> objFnEfficiency(pSol_3D_suction);
            objFnEfficiency.setEfficiencyWeight(1.0);//objFnEfficiency.setEfficiencyWeight(weightEfficiency);
            objFnEfficiency.initialize(angularVelocity, allOflowRate(nOpen), gravity, numBlades, pSol_3D_suction.parDim(), density);
            //real_t head;
            objFnEfficiency.setHead(H);
            //head = objFnEfficiency.evaluateHeadFromLoftedPressure(basisInlet, pIn_3D, velocity_relative_x_opt, gravity);
            //head -= objFnEfficiency.evaluateHeadFromLoftedPressure(basisOutlet, pOut_3D, velocity_relative_x_opt, gravity);
            //objFnEfficiency.setHead(head);
            T objFn3D_suction = objFnEfficiency.evaluateEfficiencyFromLoftedPressure(basisSuction);
            //T efficiency_suction = objFnEfficiency.getEfficiency();
            objFnEfficiency.updateSolution(pSol_3D_pressure);
            T objFn3D_pressure = objFnEfficiency.evaluateEfficiencyFromLoftedPressure(basisPressure, -1);
            //T efficiency_pressure = objFnEfficiency.getEfficiency();
            T objFn3D = objFn3D_pressure + objFn3D_suction;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (objFn3D < 0){objFn3D = 1e+10;}
            objFsw += weightEfficiency/objFn3D;
            fish.open[nOpen].f_eff = objFn3D;
        }
    }
    return objFsw;
}


private:

};
